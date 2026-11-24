#!/usr/bin/env python3
"""Sweep the scaling level and measure security / utility / compatibility.

For each setting: encrypt, decrypt as every analyst (ours + reference crypt4gh
CLI if installed), then report
  leak     aligned bases outside the analyst's regions / original such bases
  cnv_auc  AUROC of |log2(window depth / median)| for CNV vs normal 1 kb windows (0.5 = no signal)
  baf_auc  AUROC of |AF - 0.5| for het SNPs in the duplication vs elsewhere (0.5 = no signal)
  recall   bcftools SNP/indel recall outside CNVs
Usage: python3 scripts/evaluate.py [--depths raw,5,10,20,30,10q]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pysam
from sklearn.metrics import roc_auc_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from cryptbam import c4gh                                       # noqa: E402
from cryptbam.pipeline import encrypt, load_policy               # noqa: E402
from cryptbam.regions import build_atoms, readable               # noqa: E402

ANALYSTS = ["breast_analyst", "ovarian_analyst", "thal_analyst"]
sh = lambda cmd: subprocess.run(cmd, shell=True, capture_output=True, text=True)  # noqa: E731


def depth(bam, c, s, e):
    d = np.zeros(e - s + 1, int)
    for r in bam.fetch(c, s, e):
        if not (r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary):
            for a, b in r.get_blocks():
                d[max(a, s) - s] += 1
                d[min(b, e) - s] -= 1
    return np.cumsum(d[:-1])


def leaked_bases(path, mask):
    out = 0
    with pysam.AlignmentFile(path) as b:
        for r in b.fetch(until_eof=True):
            if not (r.is_unmapped or r.is_duplicate or r.is_secondary or r.is_supplementary):
                out += sum(int((~mask[r.reference_name][a:e]).sum()) for a, e in r.get_blocks())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depths", default="raw,5,10,20,30,10q")
    a = ap.parse_args()
    D, OUT = f"{ROOT}/data", f"{ROOT}/out"
    src = pysam.AlignmentFile(f"{D}/sample.bam")
    contigs = list(zip(src.references, src.lengths))
    truth = [l.rstrip("\n").split("\t") for l in open(f"{D}/truth_variants.tsv") if l[0] != "#"]
    cnvs = [l.split("\t")[:4] for l in open(f"{D}/truth_cnv.bed")]
    cnvs = [(c, int(s), int(e), k) for c, s, e, k in cnvs]
    in_cnv = lambda c, p, pad=0, kinds=("het_del", "hom_del", "dup"): any(  # noqa: E731
        cc == c and s - pad <= p < e + pad and k in kinds for cc, s, e, k in cnvs)

    pol = load_policy(f"{D}/policy.json")
    atoms = build_atoms({g: s["bed"] for g, s in pol["gene_sets"].items()}, contigs)
    classes = sorted({x.cls for x in atoms})
    regions, mask = {}, {}
    for n in ANALYSTS:
        cl = readable(pol["recipients"][n]["reads"], classes)
        regions[n] = [(x.contig, x.start, x.end) for x in atoms if x.cls in cl]
        mask[n] = {c: np.zeros(L, bool) for c, L in contigs}
        for c, s, e in regions[n]:
            mask[n][c][s:e] = True
        os.makedirs(OUT, exist_ok=True)
        with open(f"{OUT}/{n}.bed", "w") as f:
            f.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in regions[n])
    orig_out = {n: leaked_bases(f"{D}/sample.bam", mask[n]) for n in ANALYSTS}

    rows = []
    for setting in a.depths.split(","):
        k, quant = (0, False) if setting == "raw" else (int(setting.rstrip("q")), setting.endswith("q"))
        for gs in pol["gene_sets"].values():
            gs.update(depth=k, quantize_af=quant, softclip="drop")
        od = f"{OUT}/{setting}"
        os.makedirs(od, exist_ok=True)
        encrypt(f"{D}/sample.bam", pol, f"{od}/sample", log=lambda *_: None)
        enc = f"{od}/sample.bam.c4gh"
        res = {"setting": setting, "MB": round(os.path.getsize(enc) / 1e6, 2), "leak": 0,
               "ref_cli_identical": True, "samtools_ok": True, "recall": []}
        cnv_pts, baf_pts = [], []
        for n in ANALYSTS:
            bam = f"{od}/{n}.bam"
            with open(enc, "rb") as f, open(bam, "wb") as o:
                c4gh.decrypt(f, o, c4gh.load_sk(f"{D}/keys/{n}.sec"))
            if shutil.which("crypt4gh"):
                ref = sh(f"crypt4gh decrypt --sk {D}/keys/{n}.sec < {enc} | cmp -s - {bam}")
                res["ref_cli_identical"] &= ref.returncode == 0
            res["samtools_ok"] &= sh(f"samtools quickcheck {bam} && samtools index {bam}").returncode == 0
            res["leak"] = max(res["leak"], leaked_bases(bam, mask[n]) / max(orig_out[n], 1))
            b = pysam.AlignmentFile(bam)
            for c, s, e in regions[n]:                     # coverage -> CNV windows
                d = depth(b, c, s, e)
                med = np.median(d)
                for w in range(s, e - 999, 1000):
                    cnv_pts.append((any(cc == c and cs < w + 1000 and ce > w for cc, cs, ce, _ in cnvs),
                                    abs(np.log2((d[w - s:w - s + 1000].mean() + .5) / (med + .5)))))
            for c, p, ref, alt, gt, typ in truth:          # allelic imbalance at het SNPs
                p = int(p)
                if typ != "SNP" or gt not in ("1|0", "0|1") or not mask[n][c][p] \
                        or in_cnv(c, p, kinds=("het_del", "hom_del")):
                    continue
                cov = b.count_coverage(c, p, p + 1, quality_threshold=0)
                tot = sum(x[0] for x in cov)
                if tot:
                    baf_pts.append((in_cnv(c, p, kinds=("dup",)), abs(cov["ACGT".index(alt)][0] / tot - .5)))
            vcf = f"{od}/{n}.vcf.gz"                       # utility: bcftools recall
            sh(f"bcftools mpileup -f {D}/ref.fa -R {OUT}/{n}.bed {bam} -Ou | "
               f"bcftools call -mv -Oz -o {vcf} && bcftools index -f {vcf}")
            calls = {(r.chrom, r.pos - 1, r.ref, r.alts[0]) for r in pysam.VariantFile(vcf)}
            want = [(c, int(p), r, al) for c, p, r, al, *_ in truth
                    if mask[n][c][int(p)] and not in_cnv(c, int(p), 200)]
            res["recall"].append(sum(w in calls for w in want) / len(want))
        res["recall"] = round(float(np.mean(res["recall"])), 3)
        res["cnv_auc"] = round(roc_auc_score(*zip(*cnv_pts)), 3)
        res["baf_auc"] = round(roc_auc_score(*zip(*baf_pts)), 3)
        res["leak"] = f"{res['leak']:.4%}"
        rows.append(res)
        print(json.dumps(res), flush=True)

    keys = ["setting", "MB", "leak", "cnv_auc", "baf_auc", "recall", "ref_cli_identical", "samtools_ok"]
    md = ["| " + " | ".join(keys) + " |", "|" + "---|" * len(keys)]
    md += ["| " + " | ".join(str(r[k]) for k in keys) + " |" for r in rows]
    open(f"{OUT}/results.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
