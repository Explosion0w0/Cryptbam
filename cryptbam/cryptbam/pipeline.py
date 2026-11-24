"""BAM -> one region-keyed Crypt4GH BAM + one Crypt4GH bedGraph of true depth.

Plaintext layout (both files):  [_public header] [atom 1] ... [atom n] [_public EOF]
Every run = 1 junk byte + content + padding up to whole 64 KiB segments, so one
segment never mixes two keys and any subset of runs is still a valid BAM.

Policy (JSON, paths relative to the file):
{
  "gene_sets": {"breast_cancer": {"bed": "beds/breast_cancer.bed", "depth": 10},
                "a-thalassemia": {"bed": "beds/a-thalassemia.bed", "depth": 0, "softclip": "inside"}},
  "recipients": {"breast_analyst": {"pubkey": "keys/b.pub", "reads": ["breast_cancer"], "coverage": []},
                 "owner": {"pubkey": "keys/o.pub", "reads": ["*"], "coverage": ["*"], "background": true}}
}
Per gene set:  depth       0 = real reads (only clipped), K > 0 = K synthetic reads per position
               quantize_af round allele fractions to 0 / 0.5 / 1 (hides allelic imbalance)
               softclip    "drop" (default) | "inside"
A region shared by several gene sets uses the most protective setting.
"""
import gzip
import json
import os
import tempfile

import pysam

from . import bgzf, c4gh
from .regions import BACKGROUND, PUBLIC, build_atoms, read_bed, readable
from .transform import bedgraph, clip_read, consensus_reads, primary, true_depth


def load_policy(path, depth_overrides=None):
    base = os.path.dirname(os.path.abspath(path))
    pol = json.load(open(path))
    for name, gs in pol["gene_sets"].items():
        gs["bed"] = os.path.join(base, gs["bed"])
        if name in (depth_overrides or {}):
            gs["depth"] = depth_overrides[name]
    for r in pol["recipients"].values():
        r["pubkey"] = os.path.join(base, r["pubkey"])
    return pol


def scaling(gene_sets, labels, default_depth):
    """Combine the settings of every gene set covering a region (most protective wins)."""
    opts = [gene_sets[g] for g in labels]
    ks = [o.get("depth", default_depth) for o in opts]
    k = min([x for x in ks if x > 0], default=0)
    return (k, any(o.get("quantize_af") for o in opts),
            "inside" if opts and all(o.get("softclip") == "inside" for o in opts) else "drop")


def text_fill(n):
    """n bytes of bedGraph comment lines."""
    lines = []
    while n:
        k = min(n, 1000) if n - min(n, 1000) != 1 else 998
        lines.append(b"#" + b" " * (k - 2) + b"\n")
        n -= k
    return b"".join(lines)


def encrypt(bam_path, policy, out, default_depth=10, read_len=150, background=True, log=print):
    bam = pysam.AlignmentFile(bam_path)
    h = bam.header.to_dict()           # drop @PG CL (paths) and @CO (free text)
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": h["SQ"], "RG": h.get("RG", []),
         "PG": [{"ID": "cryptbam", "PN": "cryptbam"}]})
    gene_sets = policy["gene_sets"]
    beds = {g: s["bed"] for g, s in gene_sets.items()}
    # Only tile contigs that have reads or appear in a BED: a real GRCh38 header has
    # thousands of contigs and every atom costs at least one 64 KiB segment.
    used = {s.contig for s in bam.get_index_statistics() if s.mapped} | \
           {c for p in beds.values() for c, _, _ in read_bed(p)}
    atoms = build_atoms(beds, [(c, n) for c, n in zip(bam.references, bam.lengths) if c in used], log)
    tmp = os.path.join(tempfile.mkdtemp(), "x.bam")

    def raw_bam(records):              # uncompressed BAM bytes (header + records)
        with pysam.AlignmentFile(tmp, "wbu", header=header) as o:
            for r in records:
                o.write(r)
        return gzip.decompress(open(tmp, "rb").read())
    hdr = raw_bam([])

    bam_runs = [(PUBLIC, bgzf.pad_run(bgzf.compress(hdr), c4gh.SEG))]
    bdg_runs = [(PUBLIC, bgzf.pad_run(b"track type=bedGraph name=true_depth\n", c4gh.SEG, text_fill, 2))]
    manifest = []
    for i, a in enumerate(atoms):
        if a.cls == BACKGROUND and not background:
            continue
        k, quant, soft = scaling(gene_sets, a.labels, default_depth)
        if k and a.cls != BACKGROUND:
            recs = consensus_reads(bam, header, a.contig, a.start, a.end, k, read_len, quant, f"syn{i}")
        else:
            recs = [clip_read(r, a.start, a.end, header, soft) for r in bam.fetch(a.contig, a.start, a.end)
                    if primary(r)]
            recs = sorted((r for r in recs if r), key=lambda r: r.reference_start)
        bam_runs.append((a.cls, bgzf.pad_run(bgzf.compress(raw_bam(recs)[len(hdr):]), c4gh.SEG)))
        bdg, total = [], 0             # 1 Mb chunks: a whole chromosome would not fit in memory
        for s in range(a.start, a.end, 1 << 20):
            dep = true_depth(bam, a.contig, s, min(s + (1 << 20), a.end))
            bdg.append(bedgraph(a.contig, s, dep))
            total += int(dep.sum())
        bdg_runs.append((a.cls, bgzf.pad_run(b"".join(bdg), c4gh.SEG, text_fill, 2)))
        manifest.append({"region": f"{a.contig}:{a.start}-{a.end}", "class": a.cls, "depth": k,
                         "records": len(recs), "bytes": len(bam_runs[-1][1]),
                         "true_depth": round(total / (a.end - a.start), 2)})
        log(f"{manifest[-1]['region']:>28} {a.cls:<16} depth={k or 'raw'} records={len(recs)}")
    bam_runs.append((PUBLIC, b"\0" + bgzf.EOF))
    bdg_runs.append((PUBLIC, b"\0#end\n"))

    classes = sorted({a.cls for a in atoms})
    access = {}
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    for suffix, runs, field in ((".bam.c4gh", bam_runs, "reads"), (".depth.bdg.c4gh", bdg_runs, "coverage")):
        keys = {c: os.urandom(32) for c in classes + [PUBLIC]}       # fresh keys per file and class
        rcpt = []
        for name, r in policy["recipients"].items():
            cl = readable(r.get(field, []), classes, r.get("background", False))
            access.setdefault(name, {})[field] = sorted(cl)
            if cl:
                rcpt.append((c4gh.load_pk(r["pubkey"]), cl | {PUBLIC}))
        with open(out + suffix, "wb") as f:
            c4gh.write(f, runs, keys, rcpt)
    with open(out + ".owner_manifest.json", "w") as f:     # owner-only: reveals per-region record counts
        json.dump({"access": access, "atoms": manifest}, f, indent=1)
    return access, manifest
