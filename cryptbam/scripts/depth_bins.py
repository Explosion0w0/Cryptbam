#!/usr/bin/env python3
"""Mean depth per bin for BAMs or bedGraphs, plus the deletion/flank ratio.

Usage: python3 scripts/depth_bins.py REGION BIN DEL_START-DEL_END LABEL=FILE ...
  e.g. depth_bins.py chr16:140000-200000 5000 165000-185000 plain=a.bam cryptbam=b.dec.bam
BAM depth counts aligned bases like `samtools depth` (N bases included); the
"N%" column shows how much of that depth is N, i.e. "no data" written into reads.
"""
import sys

import numpy as np
import pysam


def bam_depth(path, c, s, e):
    d, n = np.zeros(e - s + 1), np.zeros(e - s + 1)
    for r in pysam.AlignmentFile(path).fetch(c, s, e):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        seq = r.query_sequence or ""
        for q, p in r.get_aligned_pairs(matches_only=True):
            if s <= p < e:
                d[p - s] += 1
                n[p - s] += seq[q] == "N"
    return d[:-1], n[:-1]


def bdg_depth(path, c, s, e):
    d = np.zeros(e - s)
    for line in open(path):
        f = line.split("\t")
        if line[0] != "#" and f[0] == c and len(f) == 4:
            a, b = max(int(f[1]), s), min(int(f[2]), e)
            if b > a:
                d[a - s:b - s] = float(f[3])
    return d, np.zeros(e - s)


c, rng = sys.argv[1].split(":")
s, e = map(int, rng.replace(",", "").split("-"))
w = int(sys.argv[2])
ds, de = map(int, sys.argv[3].split("-"))
print(f"{c}:{s}-{e}, {w // 1000} kb bins; deletion window {ds}-{de} marked with *")
print(f"{'':16}" + "".join(f"{'*' if ds <= b < de else ' ':>5}" for b in range(s, e, w)) + "   del/flank  N%")
for arg in sys.argv[4:]:
    label, path = arg.split("=", 1)
    d, n = (bdg_depth if path.endswith(".bdg") else bam_depth)(path, c, s, e)
    pos = np.arange(s, e)
    inside = (pos >= ds) & (pos < de)
    ratio = d[inside].mean() / max(d[~inside].mean(), 1e-9)
    bins = "".join(f"{d[i:i + w].mean():5.0f}" for i in range(0, e - s, w))
    npct = 100 * n[inside].sum() / max(d[inside].sum(), 1)
    flag = "  <- deletion visible" if ratio < .75 else ""
    print(f"{label:16}{bins}   {ratio:9.2f}  {npct:3.0f}{flag}")
