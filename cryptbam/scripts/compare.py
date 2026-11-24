#!/usr/bin/env python3
"""Demo: what can an attacker learn from
   (1) a plain BAM, (2) a standard single-key Crypt4GH file, (3) a cryptbam file?

Attackers:  outsider       - has the file but no key
            breast_analyst - legitimately allowed to read only the breast gene set
Run after scripts/simulate.py and keygen (see run_demo.sh).
"""
import os
import sys

import numpy as np
import pysam

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from cryptbam import c4gh                       # noqa: E402
from cryptbam.pipeline import encrypt, load_policy  # noqa: E402

D, OUT = f"{ROOT}/data", f"{ROOT}/out/compare"
BRCA1_DEL, BRCA1_FLANK = ("chr17", 43070000, 43080000), ("chr17", 43050000, 43060000)
SEA_DEL, SEA_FLANK = ("chr16", 165000, 185000), ("chr16", 140000, 160000)     # alpha-thalassemia --SEA
ALLOWED = {"chr13": (32315474, 32400266), "chr17": (43044295, 43125483)}     # breast_cancer.bed


def mean_depth(bam, c, s, e):
    return float(np.sum(bam.count_coverage(c, s, e, quality_threshold=0)) / (e - s))


def insider_view(path):
    """What the breast analyst can extract from the BAM it ends up with."""
    b = pysam.AlignmentFile(path)
    other = mates_out = 0
    for r in b.fetch(until_eof=True):
        if r.is_unmapped:
            continue
        s, e = ALLOWED.get(r.reference_name, (0, 0))
        other += sum(max(0, min(be, s) - bs) + max(0, be - max(bs, e)) for bs, be in r.get_blocks())
        if r.is_paired and r.next_reference_id >= 0:
            ms, me = ALLOWED.get(b.get_reference_name(r.next_reference_id), (0, 0))
            mates_out += not ms <= r.next_reference_start < me
    hdr = str(b.header)
    ratio = mean_depth(b, *BRCA1_DEL) / max(mean_depth(b, *BRCA1_FLANK), 1e-9)
    alpha = b.count("chr16", 140000, 200000) if "chr16" in b.references else 0
    sea = mean_depth(b, *SEA_DEL) / max(mean_depth(b, *SEA_FLANK), 1e-9) if alpha else None
    return {"unauthorized bases readable": f"{other:,}",
            "alpha-globin (chr16) reads": f"{alpha:,}",
            "alpha-thal carrier depth ratio": "no access" if sea is None else
                f"{sea:.2f}" + ("  <- carrier visible" if sea < .75 else ""),
            "mates pointing outside": f"{mates_out:,}",
            "patient note / command line": "yes" if "@CO" in hdr or "CL:" in hdr else "no",
            "BRCA1 deletion depth ratio": f"{ratio:.2f}" + ("  <- deletion visible" if ratio < .75 else "")}


def main():
    os.makedirs(OUT, exist_ok=True)
    sk = c4gh.load_sk(f"{D}/keys/breast_analyst.sec")
    pk = c4gh.load_pk(f"{D}/keys/breast_analyst.pub")
    plain = f"{D}/sample.bam"

    # (2) standard Crypt4GH: one data key for the whole file
    std = f"{OUT}/standard.bam.c4gh"
    data = open(plain, "rb").read()
    with open(std, "wb") as f:
        c4gh.write(f, [("all", b"\0" + data)], {"all": os.urandom(32)}, [(pk, {"all"})])
    with open(std, "rb") as f, open(f"{OUT}/standard.dec.bam", "wb") as o:
        c4gh.decrypt(f, o, sk)
    pysam.index(f"{OUT}/standard.dec.bam")

    # (3) cryptbam: per-gene-set keys + depth scaling
    encrypt(plain, load_policy(f"{D}/policy.json"), f"{OUT}/cryptbam", log=lambda *_: None)
    with open(f"{OUT}/cryptbam.bam.c4gh", "rb") as f, open(f"{OUT}/cryptbam.dec.bam", "wb") as o:
        st = c4gh.decrypt(f, o, sk)
    pysam.index(f"{OUT}/cryptbam.dec.bam")

    cols = {"plain BAM": plain, "standard Crypt4GH": f"{OUT}/standard.dec.bam",
            "cryptbam": f"{OUT}/cryptbam.dec.bam"}
    outsider = {"plain BAM": "everything (no key needed)",
                "standard Crypt4GH": f"file size only ({os.path.getsize(std) / 1e6:.1f} MB)",
                "cryptbam": f"file size only ({os.path.getsize(OUT + '/cryptbam.bam.c4gh') / 1e6:.1f} MB)"}
    views = {k: insider_view(v) for k, v in cols.items()}
    rows = [("[outsider] can read", outsider)] + \
           [(f"[breast analyst] {m}", {k: views[k][m] for k in cols}) for m in views["plain BAM"]]
    w = max(len(r[0]) for r in rows)
    print(f"{'':{w}}  " + "  ".join(f"{k:<30}" for k in cols))
    for name, vals in rows:
        print(f"{name:{w}}  " + "  ".join(f"{vals[k]:<30}" for k in cols))
    print(f"\ncryptbam: breast analyst decrypted {st['decrypted']} segments, "
          f"skipped {st['skipped']} it has no key for.")


if __name__ == "__main__":
    main()
