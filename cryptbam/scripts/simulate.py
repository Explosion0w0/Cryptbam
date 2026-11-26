#!/usr/bin/env python3
"""Simulate a small paired-end WGS-like BAM with *known* variants so that
leakage and utility can be measured exactly.

Contigs keep hg38 names/coordinates but are truncated just after the
simulated windows (BRCA2 / alpha-globin / BRCA1 neighbourhoods); sequence is random.
Ground truth written:  truth_variants.tsv, truth_cnv.bed, ref.fa(+.fai)
"""
import argparse
import os
import random

import numpy as np
import pysam

CONTIGS = [("chr13", 32_450_000), ("chr16", 300_000), ("chr17", 43_250_000)]
WINDOWS = {"chr13": (32_300_000, 32_420_000), "chr16": (100_000, 240_000),
           "chr17": (43_000_000, 43_200_000)}
# kind, contig, start, end, copies on (hap0, hap1)    (0-based half-open)
CNVS = [
    ("het_del", "chr17", 43_070_000, 43_080_000, (1, 0)),   # inside BRCA1  (breast+ovarian)
    ("dup",     "chr13", 32_350_000, 32_360_000, (2, 1)),   # inside BRCA2  (breast only)
    ("het_del", "chr16", 165_000,    185_000,    (1, 0)),   # --SEA alpha-thalassemia carrier (a-thalassemia)
    ("het_del", "chr17", 43_150_000, 43_160_000, (1, 0)),   # background (nobody but owner)
]
# small indels: contig, pos(0-based, variant placed *before* this ref base), ref_len, alt_ins, haps
INDELS = [
    ("chr17", 43_094_000, 3, "", (1, 0)),                      # 3 bp deletion, het
    ("chr13", 32_340_500, 0, "ACGTTG", (0, 1)),                # 6 bp insertion, het
    ("chr16", 150_000, 5, "", (1, 0)),                         # 5 bp deletion in HBZ, het (a-thalassemia)
    ("chr17", 43_125_420, 0, "TTTAAA", (1, 1)),                # hom insertion 63 bp before BRCA1 end
]
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def build_haplotypes(ref, ws, we, contig, snps, rng):
    """Return two (seq, refmap) haplotypes for window [ws, we)."""
    haps = []
    for h in (0, 1):
        seq, rmap = [], []
        p = ws
        ins_here = {i[1]: i for i in INDELS if i[0] == contig and i[4][h]}
        cnv_here = [c for c in CNVS if c[1] == contig]
        while p < we:
            if p in ins_here:
                _, _, rl, alt, _ = ins_here[p]
                for b in alt:
                    seq.append(b)
                    rmap.append(-1)
                if rl:
                    p += rl
                    continue
            cn = 1
            for kind, _, cs, ce, copies in cnv_here:
                if cs <= p < ce:
                    cn = copies[h]
                    if cn == 2 and p == cs:          # tandem duplication: emit block twice
                        for _rep in range(2):
                            for q in range(cs, ce):
                                seq.append(snps.get((q, h), ref[q - ws]))
                                rmap.append(q)
                        p = ce
                        cn = -1
                    break
            if cn == -1:
                continue
            if cn == 0:
                p += 1
                continue
            seq.append(snps.get((p, h), ref[p - ws]))
            rmap.append(p)
            p += 1
        haps.append(("".join(seq), np.array(rmap, dtype=np.int64)))
    return haps


def align(bases_rmap):
    """CIGAR for a read (in reference orientation) from its ref map; the longest
    collinear block is aligned, the rest soft-clipped (like a local aligner)."""
    blocks, cur = [], []
    prev = None
    for i, r in enumerate(bases_rmap):
        if r >= 0 and prev is not None and not (prev < r <= prev + 51):
            blocks.append(cur)
            cur = []
        cur.append(i)
        if r >= 0:
            prev = r
    blocks.append(cur)
    best = max(blocks, key=lambda b: sum(1 for i in b if bases_rmap[i] >= 0))
    # trim insertions at block edges
    while best and bases_rmap[best[0]] < 0:
        best = best[1:]
    while best and bases_rmap[best[-1]] < 0:
        best = best[:-1]
    if sum(1 for i in best if bases_rmap[i] >= 0) < 30:
        return None
    s, e = best[0], best[-1] + 1
    cig = []

    def push(op, n):
        if cig and cig[-1][0] == op:
            cig[-1] = (op, cig[-1][1] + n)
        else:
            cig.append((op, n))
    if s:
        push(4, s)
    prev = None
    for i in range(s, e):
        r = bases_rmap[i]
        if r < 0:
            push(1, 1)
        else:
            if prev is not None and r > prev + 1:
                push(2, r - prev - 1)
            push(0, 1)
            prev = r
    if e < len(bases_rmap):
        push(4, len(bases_rmap) - e)
    return int(bases_rmap[s]), [(op, int(n)) for op, n in cig]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data")
    ap.add_argument("--depth", type=float, default=30)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--skip-cnv", type=int, action="append", default=[],
                    help="index into CNVS to leave out (repeatable), e.g. 0 = no BRCA1 deletion")
    ap.add_argument("--no-ref", action="store_true", help="do not write ref.fa (many-sample runs)")
    a = ap.parse_args()
    global CNVS
    CNVS = [c for i, c in enumerate(CNVS) if i not in a.skip_cnv]
    rng = random.Random(a.seed)
    nrng = np.random.default_rng(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    L = a.read_len

    refs, snps_all, truth = {}, {}, []
    for c, ln in CONTIGS:
        ws, we = WINDOWS[c]
        refs[c] = "".join(rng.choices("ACGT", weights=[.3, .2, .2, .3], k=we - ws))
    if not a.no_ref:                         # reference FASTA (N outside windows)
        with open(os.path.join(a.outdir, "ref.fa"), "w") as f:
            for c, ln in CONTIGS:
                ws, we = WINDOWS[c]
                f.write(f">{c}\n")
                s = "N" * ws + refs[c] + "N" * (ln - we)
                for i in range(0, ln, 80):
                    f.write(s[i:i + 80] + "\n")
        pysam.faidx(os.path.join(a.outdir, "ref.fa"))

    for c, _ in CONTIGS:
        ws, we = WINDOWS[c]
        snps = {}
        for p in range(ws + 300, we - 300, 700):
            p += rng.randrange(-200, 200)
            if any(cs - 5 <= p < ce + 5 for _, cc, cs, ce, cp in CNVS if cc == c and cp == (0, 0)):
                continue
            if any(abs(p - ip) < 80 for cc, ip, *_ in INDELS if cc == c):
                continue
            refb = refs[c][p - ws]
            alt = rng.choice([b for b in "ACGT" if b != refb])
            gt = rng.choices([(1, 0), (0, 1), (1, 1)], weights=[.35, .35, .3])[0]
            for h in (0, 1):
                if gt[h]:
                    snps[(p, h)] = alt
            truth.append((c, p, refb, alt, f"{gt[0]}|{gt[1]}", "SNP"))
        snps_all[c] = snps
    for c, p, rl, alt, gt in INDELS:
        ws = WINDOWS[c][0]
        anchor = refs[c][p - 1 - ws]
        truth.append((c, p - 1, anchor + refs[c][p - ws:p - ws + rl], anchor + alt,
                      f"{gt[0]}|{gt[1]}", "INDEL"))
    with open(os.path.join(a.outdir, "truth_variants.tsv"), "w") as f:
        f.write("#contig\tpos0\tref\talt\tgt\ttype\n")
        for t in sorted(truth):
            f.write("\t".join(map(str, t)) + "\n")
    with open(os.path.join(a.outdir, "truth_cnv.bed"), "w") as f:
        for kind, c, s, e, cp in CNVS:
            f.write(f"{c}\t{s}\t{e}\t{kind}\tCN={sum(cp)}\n")

    header = {"HD": {"VN": "1.6", "SO": "unsorted"},
              "SQ": [{"SN": c, "LN": ln} for c, ln in CONTIGS],
              "RG": [{"ID": "S1", "SM": "SIM_PATIENT_001", "PL": "ILLUMINA"}],
              "PG": [{"ID": "bwa", "PN": "bwa", "CL": "bwa mem -R @RG\\tID:S1 /secret/path/patient001.fq"}],
              "CO": ["patient=SIM_PATIENT_001 referral=hereditary-cancer"]}
    tmp = os.path.join(a.outdir, "sim.unsorted.bam")
    out = pysam.AlignmentFile(tmp, "wb", header=header)
    n_pairs = 0
    for c, _ in CONTIGS:
        ws, we = WINDOWS[c]
        tid = [x[0] for x in CONTIGS].index(c)
        for h, (hseq, hmap) in enumerate(build_haplotypes(refs[c], ws, we, c, snps_all[c], rng)):
            nfrag = int(a.depth / 2 * len(hseq) / (2 * L))
            for _ in range(nfrag):
                fl = int(min(600, max(L + 20, nrng.normal(350, 35))))
                s = rng.randrange(0, len(hseq) - fl)
                n_pairs += 1
                name = f"SIM:{c}:{n_pairs}"
                pair = []
                for side in (0, 1):
                    a0 = s if side == 0 else s + fl - L
                    bases = list(hseq[a0:a0 + L])
                    qual = [rng.randint(28, 40) for _ in range(L)]
                    for i in range(L):
                        if rng.random() < 0.002:
                            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
                            qual[i] = rng.randint(8, 20)
                    al = align(list(hmap[a0:a0 + L]))
                    pair.append((bases, qual, al, side == 1))
                dup = rng.random() < 0.02
                for copy in range(2 if dup else 1):
                    recs = []
                    first = rng.random() < 0.5
                    for side, (bases, qual, al, rev) in enumerate(pair):
                        r = pysam.AlignedSegment(out.header)
                        r.query_name = name + ("" if copy == 0 else ":dup")
                        r.query_sequence = "".join(bases)
                        r.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in qual))
                        r.flag = 0x1 | (0x40 if (side == 0) == first else 0x80)
                        if rev:
                            r.flag |= 0x10
                        if copy:
                            r.flag |= 0x400
                        if al is None:
                            r.flag |= 0x4
                        else:
                            r.reference_id, r.reference_start = tid, al[0]
                            r.cigartuples = al[1]
                            r.mapping_quality = 60
                            # NM
                            nm, rp, qp = 0, al[0], 0
                            for op, ln in al[1]:
                                if op == 0:
                                    for k in range(ln):
                                        if bases[qp + k] != refs[c][rp + k - ws]:
                                            nm += 1
                                    rp += ln; qp += ln
                                elif op == 1:
                                    nm += ln; qp += ln
                                elif op == 2:
                                    nm += ln; rp += ln
                                elif op == 4:
                                    qp += ln
                            r.set_tag("NM", int(nm))
                        r.set_tag("RG", "S1")
                        recs.append(r)
                    r1, r2 = recs
                    for x, y in ((r1, r2), (r2, r1)):
                        if y.is_unmapped and not x.is_unmapped:
                            y.reference_id, y.reference_start = x.reference_id, x.reference_start
                    for x, y in ((r1, r2), (r2, r1)):
                        x.next_reference_id, x.next_reference_start = y.reference_id, y.reference_start
                        if y.is_unmapped:
                            x.flag |= 0x8
                        if y.is_reverse:
                            x.flag |= 0x20
                    if not r1.is_unmapped and not r2.is_unmapped:
                        lo = min(r1.reference_start, r2.reference_start)
                        hi = max(r1.reference_end, r2.reference_end)
                        for x in (r1, r2):
                            x.flag |= 0x2
                            x.template_length = (hi - lo) * (1 if x.reference_start == lo else -1)
                            x.set_tag("MC", (r2 if x is r1 else r1).cigarstring)
                    out.write(r1)
                    out.write(r2)
    out.close()
    final = os.path.join(a.outdir, "sample.bam")
    pysam.sort("-o", final, tmp)
    pysam.index(final)
    os.remove(tmp)
    print(f"wrote {final}: {n_pairs} pairs")


if __name__ == "__main__":
    main()
