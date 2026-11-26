#!/usr/bin/env python3
"""Quantify cryptbam against a plain BAM and standard (single-key) Crypt4GH.

Attackers: an OUTSIDER (has the file, no key) and each INSIDER analyst (holds a
valid key for its own gene set). With standard Crypt4GH an insider decrypts the
whole original BAM; with cryptbam it gets only its own regions.

  A  simulated sample with known CNVs/variants (data/)  -> security + utility
  B  N simulated samples, half with the BRCA1 deletion  -> size-inference AUROC
  C  real data: HG00684 (alpha-thal carrier) + NA12878  -> security, utility, timing

Usage: python3 scripts/benchmark.py [--samples 8] [--skip-size] [--skip-real]
Writes out/bench/benchmark.md (+ .json).
"""
import argparse
import json
import os
import struct
import subprocess
import sys
import time

import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from cryptbam import c4gh                                     # noqa: E402
from cryptbam.pipeline import encrypt, load_policy             # noqa: E402
from cryptbam.regions import build_atoms, readable             # noqa: E402

OUT = f"{ROOT}/out/bench"
ANALYSTS = ["breast_analyst", "ovarian_analyst", "thal_analyst"]
SCHEMES = ["plain", "standard", "cryptbam"]
sh = lambda cmd: subprocess.run(cmd, shell=True, capture_output=True, text=True)  # noqa: E731
quiet = lambda *_: None                                                          # noqa: E731


# ---- helpers ---------------------------------------------------------------------
def std_encrypt(bam, pk, out):
    with open(out, "wb") as f:
        c4gh.write(f, [("all", b"\0" + open(bam, "rb").read())], {"all": os.urandom(32)}, [(c4gh.load_pk(pk), {"all"})])


def dec(path, sk, out):
    """Decrypt with sk; None if the key opens nothing."""
    try:
        with open(path, "rb") as f, open(out, "wb") as o:
            c4gh.decrypt(f, o, c4gh.load_sk(sk))
    except PermissionError:
        return None
    if out.endswith(".bam"):
        pysam.index(out)
    return out


def usable(r):
    return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate)


def spans(bam):
    s = {}
    for r in pysam.AlignmentFile(bam).fetch(until_eof=True):
        if usable(r):
            a, b = s.get(r.reference_name, (r.reference_start, r.reference_end))
            s[r.reference_name] = (min(a, r.reference_start), max(b, r.reference_end))
    return s


def depth(bam, sp):
    """{contig: per-base depth over its span}; all zeros for a view with no reads."""
    d = {c: np.zeros(e - s + 1, np.int32) for c, (s, e) in sp.items()}
    if bam:
        for r in pysam.AlignmentFile(bam).fetch(until_eof=True):
            if usable(r) and r.reference_name in sp:
                s0, e0 = sp[r.reference_name]
                for a, b in r.get_blocks():
                    a, b = max(a, s0), min(b, e0)
                    if b > a:
                        d[r.reference_name][a - s0] += 1
                        d[r.reference_name][b - s0] -= 1
    return {c: np.cumsum(x[:-1]) for c, x in d.items()}


def mask(sp, intervals):
    m = {c: np.zeros(e - s, bool) for c, (s, e) in sp.items()}
    for c, a, b in intervals:
        if c in sp:
            s0, e0 = sp[c]
            m[c][max(a, s0) - s0:max(min(b, e0), s0) - s0] = True
    return m


def bdg_depth(path, sp):
    d = {c: np.zeros(e - s, np.int32) for c, (s, e) in sp.items()}
    for line in open(path) if path else []:
        f = line.split("\t")
        if line[0] != "#" and len(f) == 4 and f[0] in sp:
            s0, e0 = sp[f[0]]
            a, b = max(int(f[1]), s0), min(int(f[2]), e0)
            if b > a:
                d[f[0]][a - s0:b - s0] = int(f[3])
    return d


def ratio(dep, sp, c, ds, de, fs, fe):
    s0 = sp[c][0]
    flank = np.r_[dep[c][fs - s0:ds - s0], dep[c][de - s0:fe - s0]].mean()
    return dep[c][ds - s0:de - s0].mean() / flank if flank else "no access"


def rnd(x):
    return round(x, 2) if isinstance(x, float) else x


def auc(y, x):
    a = roc_auc_score(y, x) if len(set(y)) == 2 else float("nan")
    return round(max(a, 1 - a), 3)       # the attacker may pick either direction


def pct(a, b):
    return f"{100 * a / b:.2f}%" if b else "n/a"


def header_packets(path):
    return struct.unpack("<I", open(path, "rb").read(16)[12:16])[0]


def policy_atoms(pol, bam):
    b = pysam.AlignmentFile(bam)
    atoms = build_atoms({g: s["bed"] for g, s in pol["gene_sets"].items()}, list(zip(b.references, b.lengths)), quiet)
    classes = sorted({a.cls for a in atoms})
    allowed = {n: [(a.contig, a.start, a.end) for a in atoms
                   if a.cls in readable(pol["recipients"][n].get("reads", []), classes)] for n in ANALYSTS}
    targets = [(a.contig, a.start, a.end) for a in atoms if a.labels]
    return allowed, targets


def insider_rows(orig, views, bdgs, allowed, sp, odep):
    """Per analyst x scheme: sequence / coverage / mate / metadata leakage + retention."""
    rows = []
    for n in ANALYSTS:
        m = mask(sp, allowed[n])
        o_out = sum(int(odep[c][~m[c]].sum()) for c in sp)
        o_cov = sum(int((odep[c][~m[c]] > 0).sum()) for c in sp)
        o_in = sum(int((odep[c][m[c]] > 0).sum()) for c in sp)
        for sch in SCHEMES:
            v = views[sch][n]
            vd = depth(v, sp)
            bd = bdg_depth(bdgs[sch].get(n), sp) if sch == "cryptbam" else vd
            b = pysam.AlignmentFile(v)
            mates = sum(1 for r in b.fetch(until_eof=True) if usable(r) and r.is_paired and not r.mate_is_unmapped
                        and not any(c == r.next_reference_name and s <= r.next_reference_start < e
                                    for c, s, e in allowed[n]))
            hdr = str(b.header)
            rows.append({"analyst": n, "scheme": sch,
                         "seq_leak": pct(sum(int(vd[c][~m[c]].sum()) for c in sp), o_out),
                         "cov_leak": pct(sum(int(((vd[c][~m[c]] > 0) | (bd[c][~m[c]] > 0)).sum()) for c in sp), o_cov),
                         "mate_leak": mates,
                         "metadata": "yes" if ("CL:" in hdr or "@CO" in hdr) else "no",
                         "auth_positions_kept": pct(sum(int((vd[c][m[c]] > 0).sum()) for c in sp), o_in)})
    return rows


def calls(bam, ref, bed):
    r = sh(f"bcftools mpileup -f {ref} -R {bed} {bam} 2>/dev/null | bcftools call -mv 2>/dev/null | grep -v '^#'")
    return {(f[0], int(f[1]) - 1, f[3], f[4]) for f in (l.split("\t") for l in r.stdout.splitlines())}


def bed(intervals, path):
    with open(path, "w") as f:
        f.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in intervals)
    return path


def build_views(bam, pol, keys, od, std_for="breast_analyst"):
    """Encrypt both ways, decrypt as every analyst; record sizes and timings."""
    os.makedirs(od, exist_ok=True)
    perf = {"input_MB": os.path.getsize(bam) / 1e6}
    t = time.time(); std_encrypt(bam, f"{keys}/{std_for}.pub", f"{od}/standard.c4gh"); perf["std_encrypt_s"] = time.time() - t
    t = time.time(); dec(f"{od}/standard.c4gh", f"{keys}/{std_for}.sec", f"{od}/standard.dec.bam"); perf["std_decrypt_s"] = time.time() - t
    t = time.time(); encrypt(bam, pol, f"{od}/cb", log=quiet); perf["cb_encrypt_s"] = time.time() - t
    views = {"plain": {}, "standard": {}, "cryptbam": {}}
    bdgs = {"cryptbam": {}}
    for n in ANALYSTS:
        views["plain"][n] = bam
        views["standard"][n] = f"{od}/standard.dec.bam"   # every insider gets the whole file
        t = time.time(); views["cryptbam"][n] = dec(f"{od}/cb.bam.c4gh", f"{keys}/{n}.sec", f"{od}/{n}.bam")
        perf[f"cb_decrypt_{n}_s"] = time.time() - t
        bdgs["cryptbam"][n] = dec(f"{od}/cb.depth.bdg.c4gh", f"{keys}/{n}.sec", f"{od}/{n}.bdg")
    perf.update(std_MB=os.path.getsize(f"{od}/standard.c4gh") / 1e6,
                cb_reads_MB=os.path.getsize(f"{od}/cb.bam.c4gh") / 1e6,
                cb_depth_MB=os.path.getsize(f"{od}/cb.depth.bdg.c4gh") / 1e6,
                std_packets=header_packets(f"{od}/standard.c4gh"), cb_packets=header_packets(f"{od}/cb.bam.c4gh"),
                outsider_std=dec(f"{od}/standard.c4gh", f"{keys}/outsider.sec", f"{od}/x") is None,
                outsider_cb=dec(f"{od}/cb.bam.c4gh", f"{keys}/outsider.sec", f"{od}/x") is None)
    return views, bdgs, perf


# ---- A: simulated sample with ground truth ---------------------------------------
def part_a(R):
    D = f"{ROOT}/data"
    pol = load_policy(f"{D}/policy.json")
    bam = f"{D}/sample.bam"
    for n in ANALYSTS + ["outsider"]:
        if not os.path.exists(f"{D}/keys/{n}.sec"):
            c4gh.keygen(f"{D}/keys/{n}.sec", f"{D}/keys/{n}.pub")
    views, bdgs, perf = build_views(bam, pol, f"{D}/keys", f"{OUT}/sim")
    allowed, targets = policy_atoms(pol, bam)
    sp = spans(bam)
    odep = depth(bam, sp)
    R["A_insider"] = insider_rows(bam, views, bdgs, allowed, sp, odep)
    R["A_perf"] = perf

    cnvs = [(c, int(s), int(e), k) for c, s, e, k, *_ in (l.split("\t") for l in open(f"{D}/truth_cnv.bed"))]
    truth = [l.rstrip("\n").split("\t") for l in open(f"{D}/truth_variants.tsv") if l[0] != "#"]
    in_cnv = lambda c, p, pad=0: any(cc == c and s - pad <= p < e + pad for cc, s, e, _ in cnvs)  # noqa: E731
    rows = []
    for n in ANALYSTS:
        bedf = bed(allowed[n], f"{OUT}/sim/{n}.bed")
        want = {(c, int(p), r, a) for c, p, r, a, *_ in truth
                if any(cc == c and s <= int(p) < e for cc, s, e in allowed[n]) and not in_cnv(c, int(p), 200)}
        for sch in SCHEMES:
            v = views[sch][n]
            vd = depth(v, sp)
            cnv_pts = []                     # 1 kb windows over ALL gene-set regions, scored from this view
            for c, s, e in targets:
                s0 = sp[c][0]
                d = vd[c][s - s0:e - s0]
                med = np.median(d)
                for w in range(0, e - s - 999, 1000):
                    cnv_pts.append((any(cc == c and cs < s + w + 1000 and ce > s + w for cc, cs, ce, _ in cnvs),
                                    abs(np.log2((d[w:w + 1000].mean() + .5) / (med + .5)))))
            baf_pts = []
            b = pysam.AlignmentFile(v)
            for c, p, ref, alt, gt, typ in truth:
                p = int(p)
                if typ == "SNP" and gt in ("1|0", "0|1") and any(cc == c and s <= p < e for cc, s, e in targets) \
                        and not any(cc == c and s <= p < e and k != "dup" for cc, s, e, k in cnvs):
                    cov = b.count_coverage(c, p, p + 1, quality_threshold=0)
                    tot = sum(x[0] for x in cov)
                    baf_pts.append((any(cc == c and s <= p < e and k == "dup" for cc, s, e, k in cnvs),
                                    abs(cov["ACGT".index(alt)][0] / tot - .5) if tot else 0))
            got = calls(v, f"{D}/ref.fa", bedf)
            rows.append({"analyst": n, "scheme": sch,
                         "cnv_auroc_depth": auc(*zip(*cnv_pts)), "cnv_auroc_allele": auc(*zip(*baf_pts)),
                         "brca1_del_ratio": rnd(ratio(vd, sp, "chr17", 43070000, 43080000, 43045000, 43125000)),
                         "sea_del_ratio": rnd(ratio(vd, sp, "chr16", 165000, 185000, 140000, 200000)),
                         "recall": pct(len(want & got), len(want))})
    R["A_inference"] = rows


# ---- B: size inference over many samples ------------------------------------------
def part_b(R, n):
    D = f"{ROOT}/data"
    keys = f"{D}/keys"
    feats = {k: [] for k in ("std_outsider_size", "cb10_outsider_size", "cb10_thal_skipped",
                             "cb0_outsider_size", "cb0_thal_skipped", "plain_insider_ratio", "cb10_breast_ratio")}
    y = []
    for i in range(n):
        od = f"{OUT}/size/s{i}"
        has = i % 2 == 0
        if not os.path.exists(f"{od}/sample.bam.bai"):
            sh(f"python3 {ROOT}/scripts/simulate.py --outdir {od} --seed {100 + i} --no-ref" + ("" if has else " --skip-cnv 0"))
        bam = f"{od}/sample.bam"
        y.append(int(has))
        std_encrypt(bam, f"{keys}/breast_analyst.pub", f"{od}/std.c4gh")
        feats["std_outsider_size"].append(os.path.getsize(f"{od}/std.c4gh"))
        for k in (10, 0):
            encrypt(bam, load_policy(f"{D}/policy.json", {g: k for g in ("breast_cancer", "ovarian_cancer", "a-thalassemia")}),
                    f"{od}/cb{k}", log=quiet)
            feats[f"cb{k}_outsider_size"].append(os.path.getsize(f"{od}/cb{k}.bam.c4gh"))
            with open(f"{od}/cb{k}.bam.c4gh", "rb") as f:          # what the thal analyst's edit list reveals
                _, edits = c4gh.read_header(f, c4gh.load_sk(f"{keys}/thal_analyst.sec"))
            feats[f"cb{k}_thal_skipped"].append(sum(edits[0::2]))
        sp = spans(bam)
        feats["plain_insider_ratio"].append(ratio(depth(bam, sp), sp, "chr17", 43070000, 43080000, 43045000, 43125000))
        v = dec(f"{od}/cb10.bam.c4gh", f"{keys}/breast_analyst.sec", f"{od}/breast10.bam")
        feats["cb10_breast_ratio"].append(ratio(depth(v, sp), sp, "chr17", 43070000, 43080000, 43045000, 43125000))
    R["B_size"] = {"samples": n, "with_deletion": sum(y), **{k: auc(y, v) for k, v in feats.items()}}


# ---- C: real data -----------------------------------------------------------------
def part_c(R):
    demo = os.path.join(os.path.dirname(ROOT), "cryptbam_demo")
    pol = load_policy(f"{demo}/policy.json")
    bam = f"{demo}/seq/breast_cancer_patient_seq.bam"
    views, bdgs, perf = build_views(bam, pol, f"{demo}/keys", f"{OUT}/real")
    allowed, targets = policy_atoms(pol, bam)
    sp = spans(bam)
    R["C_insider"] = insider_rows(bam, views, bdgs, allowed, sp, depth(bam, sp))
    R["C_perf"] = perf
    # utility: concordance with the calls made on the plain BAM, per analyst region
    util = []
    ref = f"{demo}/seq/ref_genome.fa.gz"
    for n in ANALYSTS:
        bedf = bed(allowed[n], f"{OUT}/real/{n}.bed")
        base = calls(bam, ref, bedf)
        got = calls(views["cryptbam"][n], ref, bedf)
        util.append({"analyst": n, "plain_calls": len(base), "cryptbam_calls": len(got), "shared": len(base & got)})
    R["C_utility"] = util
    # alpha-thalassemia carrier: HG00684 vs normal NA12878, as seen by the thal analyst
    na = f"{ROOT}/demo/data/NA12878.genes.bam"
    car = []
    for name, b in (("HG00684 (carrier)", bam), ("NA12878 (normal)", na)):
        if not os.path.exists(b):
            continue
        od = f"{OUT}/real/{name.split()[0]}"
        os.makedirs(od, exist_ok=True)
        encrypt(b, pol, f"{od}/cb", log=quiet)
        s = {"chr16": (130000, 210000)}
        thal = dec(f"{od}/cb.bam.c4gh", f"{demo}/keys/thal_analyst.sec", f"{od}/thal.bam")
        cnv = dec(f"{od}/cb.depth.bdg.c4gh", f"{demo}/keys/cnv_analyst.sec", f"{od}/cnv.bdg")
        r = lambda d: rnd(ratio(d, s, "chr16", 165000, 185000, 140000, 200000))  # noqa: E731
        car.append({"sample": name, "plain": r(depth(b, s)), "standard": r(depth(b, s)),
                    "cryptbam_thal": r(depth(thal, s)), "cryptbam_cnv_analyst": r(bdg_depth(cnv, s))})
    R["C_carrier"] = car


# ---- report -----------------------------------------------------------------------
def table(rows, keys):
    out = ["| " + " | ".join(keys) + " |", "|" + "---|" * len(keys)]
    return out + ["| " + " | ".join(str(r[k]) for k in keys) + " |" for r in rows]


def report(R):
    md = ["# cryptbam benchmark", ""]
    for part, title in (("A", "Simulated 30x sample (known CNVs / variants)"), ("C", "Real data: HG00684 (1000 Genomes)")):
        if f"{part}_insider" not in R:
            continue
        md += [f"## {title}", "", "### Insider leakage (analyst holds a valid key)", ""]
        md += table(R[f"{part}_insider"], ["analyst", "scheme", "seq_leak", "cov_leak", "mate_leak", "metadata",
                                           "auth_positions_kept"]) + [""]
        if part == "A":
            md += ["### Inference and utility", ""]
            md += table(R["A_inference"], ["analyst", "scheme", "cnv_auroc_depth", "cnv_auroc_allele",
                                           "brca1_del_ratio", "sea_del_ratio", "recall"]) + [""]
        else:
            md += ["### Utility (calls shared with the plain BAM)", ""] + table(
                R["C_utility"], ["analyst", "plain_calls", "cryptbam_calls", "shared"]) + [""]
            md += ["### Alpha-thalassemia carrier depth ratio (deletion / flanks)", ""] + table(
                R["C_carrier"], ["sample", "plain", "standard", "cryptbam_thal", "cryptbam_cnv_analyst"]) + [""]
        p = R[f"{part}_perf"]
        md += ["### Outsider view and cost", "", "| metric | standard Crypt4GH | cryptbam |", "|---|---|---|",
               f"| outsider can decrypt | {'no' if p['outsider_std'] else 'YES'} | {'no' if p['outsider_cb'] else 'YES'} |",
               f"| header packets visible | {p['std_packets']} | {p['cb_packets']} |",
               f"| encrypted size MB (input {p['input_MB']:.2f}) | {p['std_MB']:.2f} | "
               f"{p['cb_reads_MB']:.2f} reads + {p['cb_depth_MB']:.2f} depth |",
               f"| encrypt s | {p['std_encrypt_s']:.2f} | {p['cb_encrypt_s']:.2f} |",
               f"| decrypt s (breast analyst) | {p['std_decrypt_s']:.2f} | {p['cb_decrypt_breast_analyst_s']:.2f} |", ""]
    if "B_size" in R:
        b = R["B_size"]
        md += ["## Size inference (BRCA1 deletion present / absent)", "",
               f"{b['samples']} simulated samples, {b['with_deletion']} with the deletion; AUROC = max(auc, 1-auc).", "",
               "| attacker feature | AUROC |", "|---|---|"]
        md += [f"| {k} | {v} |" for k, v in b.items() if k not in ("samples", "with_deletion")] + [""]
    open(f"{OUT}/benchmark.md", "w").write("\n".join(md))
    json.dump(R, open(f"{OUT}/benchmark.json", "w"), indent=1, default=str)
    print("\n".join(md))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=8)
    ap.add_argument("--skip-size", action="store_true")
    ap.add_argument("--skip-real", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    R = {}
    part_a(R)
    if not a.skip_real:
        part_c(R)
    if not a.skip_size:
        part_b(R, a.samples)
    report(R)


if __name__ == "__main__":
    main()
