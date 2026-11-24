#!/usr/bin/env python3
"""Demo step 6: what each identity sees after decrypting the same .bam.c4gh.

Usage: python3 scripts/identities.py DEMO_DIR    (after scripts/demo_real.sh step 4)
Prints reads / mean depth per gene, and a fingerprint of BRCA1 so you can see that
the breast and ovarian analysts get the very same bytes for the shared gene.
"""
import hashlib
import os
import sys

import pysam

D = sys.argv[1] if len(sys.argv) > 1 else "demo"
GENES = {}
for f in sorted(os.listdir(f"{D}/beds")):
    for line in open(f"{D}/beds/{f}"):
        c, s, e, name = line.split()[:4]
        GENES[name] = (c, int(s), int(e))
VIEWS = {"plain BAM": f"{D}/data/NA12878.genes.bam",
         **{n: f"{D}/out/{n}.bam" for n in ("breast_analyst", "ovarian_analyst", "thal_analyst", "owner")}}


def stats(path, c, s, e):
    b = pysam.AlignmentFile(path)
    n = b.count(c, s, e)
    dep = sum(sum(x) for x in b.count_coverage(c, s, e, quality_threshold=0)) / (e - s) if n else 0
    h = hashlib.md5("".join(r.to_string() for r in b.fetch(c, s, e)).encode()).hexdigest()[:8] if n else "-"
    return n, dep, h


rows = {v: {g: stats(p, *GENES[g]) for g in GENES} for v, p in VIEWS.items()}
w = 16
print(f"{'reads / mean depth':<{w}}" + "".join(f"{g:>18}" for g in GENES))
for v in VIEWS:
    print(f"{v:<{w}}" + "".join(
        f"{(f'{n:,} / {d:.0f}x' if n else 'no access'):>18}" for n, d, _ in rows[v].values()))
print("\nBRCA1 fingerprint (md5 of the records):")
for v in VIEWS:
    print(f"  {v:<{w}} {rows[v]['BRCA1'][2]}")
print("breast_analyst and ovarian_analyst share one key for BRCA1 -> identical records.")
