"""Read transformations applied before encryption.

clip_read        cut an alignment to [ws, we) and strip everything that can
                 leak outside it (soft clips, mate info, tags)
consensus_reads  re-emit the pileup as exactly `depth` synthetic reads per
                 position: the user-chosen scaling level (0 = keep real reads)
true_depth       real per-base depth, stored separately in the bedGraph
"""
from array import array
from collections import Counter, defaultdict

import numpy as np
import pysam

BASES = "ACGT"
LUT = np.full(256, -1, np.int8)
for i, b in enumerate(BASES):
    LUT[ord(b)] = i
DEL, NODATA = 4, 5


def primary(r):
    return not (r.is_unmapped or r.is_secondary or r.is_supplementary) and r.query_sequence


# ---- boundary read protection -------------------------------------------------
def clip_read(r, ws, we, header, softclip="drop"):
    """Part of `r` aligned inside [ws, we) as a new unpaired record, or None.
    softclip="inside" keeps a soft clip only if that read end was not cut and
    the clip projected onto the reference stays inside the window."""
    seq = r.query_sequence
    qual = r.query_qualities if r.query_qualities is not None else array("B", [30] * len(seq))
    parts, rpos, q = [], r.reference_start, 0           # parts: [op, len, query_start, ref_pos]
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            a, b = max(rpos, ws), min(rpos + n, we)
            if b > a:
                parts.append([0, b - a, q + a - rpos, a])
            rpos, q = rpos + n, q + n
        elif op == 1:
            if ws < rpos < we:
                parts.append([1, n, q, rpos])
            q += n
        elif op in (2, 3):
            a, b = max(rpos, ws), min(rpos + n, we)
            if b > a:
                parts.append([2, b - a, q, a])
            rpos += n
        elif op == 4:
            q += n
    while parts and parts[0][0]:
        parts.pop(0)
    while parts and parts[-1][0]:
        parts.pop()
    if not parts:
        return None
    if softclip == "inside":
        cig = [c for c in r.cigartuples if c[0] != 5]
        if cig[0][0] == 4 and parts[0][3] == r.reference_start and r.reference_start - cig[0][1] >= ws:
            parts.insert(0, [4, cig[0][1], 0, None])
        if cig[-1][0] == 4 and parts[-1][3] + parts[-1][1] == r.reference_end \
                and r.reference_end + cig[-1][1] <= we:
            parts.append([4, cig[-1][1], len(seq) - cig[-1][1], None])
    out = pysam.AlignedSegment(header)
    out.query_name = r.query_name + ("/1" if r.is_read1 else "/2" if r.is_read2 else "")
    out.flag = r.flag & 0x610                            # keep reverse / qcfail / duplicate only
    out.reference_id = header.get_tid(r.reference_name)
    out.reference_start = next(p[3] for p in parts if p[0] == 0)
    out.mapping_quality = r.mapping_quality
    out.cigartuples = [(p[0], p[1]) for p in parts]
    take = [(p[2], p[2] + p[1]) for p in parts if p[0] != 2]
    out.query_sequence = "".join(seq[a:b] for a, b in take)
    out.query_qualities = array("B", [x for a, b in take for x in qual[a:b]])
    if r.has_tag("RG"):
        out.set_tag("RG", r.get_tag("RG"))               # MD/NM/MC/SA/XA etc. are dropped
    return out


# ---- true depth -------------------------------------------------------------------
def true_depth(bam, contig, ws, we):
    d = np.zeros(we - ws + 1, np.int64)
    for r in bam.fetch(contig, ws, we):
        if primary(r) and not (r.is_duplicate or r.is_qcfail):
            for a, b in r.get_blocks():
                a, b = max(a, ws), min(b, we)
                if b > a:
                    d[a - ws] += 1
                    d[b - ws] -= 1
    return np.cumsum(d[:-1])


def bedgraph(contig, ws, depth):
    if not len(depth):
        return b""
    cut = np.r_[0, np.flatnonzero(np.diff(depth)) + 1, len(depth)]
    return "".join(f"{contig}\t{ws + s}\t{ws + e}\t{depth[s]}\n" for s, e in zip(cut, cut[1:])).encode()


# ---- pileup compression -----------------------------------------------------------
def pileup(bam, contig, ws, we, min_mapq=1, min_bq=13):
    """counts (n,5: A C G T del), base-quality sums (n,4), MAPQ sums (n), insertions {pos: Counter}."""
    n = we - ws
    cnt, bq, mq, ins = np.zeros((n, 5), np.int32), np.zeros((n, 4)), np.zeros(n), defaultdict(Counter)
    for r in bam.fetch(contig, ws, we):
        if not primary(r) or r.is_duplicate or r.is_qcfail or r.mapping_quality < min_mapq:
            continue
        seq = np.frombuffer(r.query_sequence.encode(), np.uint8)
        qual = np.asarray(r.query_qualities if r.query_qualities is not None else [30] * len(seq))
        rpos, q = r.reference_start, 0
        for op, ln in r.cigartuples:
            a, b = max(rpos, ws), min(rpos + ln, we)
            if op in (0, 7, 8):
                if b > a:
                    s = q + a - rpos
                    base, bqs = LUT[seq[s:s + b - a]].astype(int), qual[s:s + b - a]
                    ok = (base >= 0) & (bqs >= min_bq)
                    pos = np.arange(a - ws, b - ws)[ok]
                    np.add.at(cnt, (pos, base[ok]), 1)
                    np.add.at(bq, (pos, base[ok]), bqs[ok])
                    mq[a - ws:b - ws] += r.mapping_quality
                rpos, q = rpos + ln, q + ln
            elif op == 1:
                if ws < rpos < we:
                    ins[rpos - ws][r.query_sequence[q:q + ln]] += 1
                q += ln
            elif op in (2, 3):
                if op == 2 and b > a:
                    cnt[a - ws:b - ws, DEL] += 1
                rpos += ln
            elif op == 4:
                q += ln
    return cnt, bq, mq, ins


def layers(cnt, K, quantize):
    """(n,5) counts -> (n,K) allele per layer, largest-remainder rounding to 1/K.
    quantize: snap the alt fraction to 0 or 1/2 (hides allelic imbalance)."""
    order = np.argsort(-cnt, axis=1, kind="stable")
    c = np.take_along_axis(cnt, order, axis=1).astype(float)
    if quantize:
        het = c[:, 1] >= 0.2 * (c[:, 0] + c[:, 1])
        c = np.zeros_like(c)
        c[:, 0], c[het, 1] = 1, 1
    share = c * K / np.maximum(c.sum(1, keepdims=True), 1)
    n_al = np.floor(share).astype(int)
    rank = np.argsort(np.argsort(-(share - n_al), axis=1, kind="stable"), axis=1)
    n_al += rank < (K - n_al.sum(1))[:, None]
    cum = np.cumsum(n_al, axis=1)
    out = np.stack([np.take_along_axis(order, np.minimum((cum <= j).sum(1), 4)[:, None], 1)[:, 0]
                    for j in range(K)], axis=1)
    out[cnt.sum(1) == 0] = NODATA
    return out


def consensus_reads(bam, header, contig, ws, we, K, read_len=150, quantize=False, name="syn"):
    """Exactly K reads cover every position of [ws, we): depth carries no information.
    Positions without data get N (BQ 2, MAPQ 0); allele fractions are rounded to 1/K."""
    cnt, bq, mq, ins = pileup(bam, contig, ws, we)
    n, d = we - ws, cnt.sum(1)
    lay = layers(cnt, K, quantize)
    meanq = np.clip(np.rint(bq / np.maximum(cnt[:, :4], 1)), 2, 40).astype(int)
    meanmq = mq / np.maximum(d, 1)
    ins_at = [dict() for _ in range(K)]          # alt insertions go to the last layers
    for p, ctr in ins.items():
        f = sum(ctr.values()) / max(d[p - 1], 1)
        f = (0 if f < .2 else .5 if f < .8 else 1) if quantize else min(f, 1)
        free = [j for j in reversed(range(K)) if DEL not in (lay[p - 1, j], lay[p, j])]
        for j in free[:round(K * f)]:
            ins_at[j][p] = ctr.most_common(1)[0][0]
    reads = []
    for j in range(K):
        al, il = lay[:, j], ins_at[j]
        cuts, c = [0], (j * read_len) // K or read_len    # staggered starts per layer
        while c < n:                                       # avoid cutting at a deletion/insertion
            c0 = c
            while c < n and c - c0 < read_len and (DEL in (al[c], al[c - 1]) or c in il):
                c += 1
            cuts.append(min(c, n))
            c += read_len
        cuts.append(n)
        for k, (s, e) in enumerate(zip(cuts, cuts[1:])):
            while s < e and al[s] == DEL:
                s += 1
            while e > s and al[e - 1] == DEL:
                e -= 1
            if s >= e:
                continue
            cig, seq, qual = [], [], []
            for p in range(s, e):
                if p > s and p in il:
                    cig.append((1, len(il[p])))
                    seq.append(il[p])
                    qual += [30] * len(il[p])
                if al[p] == DEL:
                    cig.append((2, 1))
                else:
                    cig.append((0, 1))
                    seq.append("N" if al[p] == NODATA else BASES[al[p]])
                    qual.append(2 if al[p] == NODATA else int(meanq[p, al[p]]))
            has = d[s:e] > 0
            r = pysam.AlignedSegment(header)
            r.query_name, r.flag = f"{name}.{j}.{k}", 16 * (j % 2)
            r.reference_id, r.reference_start = header.get_tid(contig), ws + s
            r.mapping_quality = int(min(60, round(meanmq[s:e][has].mean()))) if has.any() else 0
            r.cigartuples = merge(cig)
            r.query_sequence, r.query_qualities = "".join(seq), array("B", qual)
            reads.append(r)
    return sorted(reads, key=lambda x: x.reference_start)


def merge(cig):
    out = []
    for op, n in cig:
        if out and out[-1][0] == op:
            out[-1] = (op, out[-1][1] + n)
        else:
            out.append((op, n))
    return out
