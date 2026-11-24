"""Gene-set BEDs -> atoms -> access classes.

An atom is a maximal interval covered by a constant set of gene sets; atoms
with the same set share one key. BRCA1 in breast.bed and ovarian.bed becomes
class "breast+ovarian": stored once, key given to both analysts.
"""
from dataclasses import dataclass

BACKGROUND, PUBLIC = "_background", "_public"


@dataclass(frozen=True)
class Atom:
    contig: str
    start: int
    end: int
    labels: frozenset

    @property
    def cls(self):
        return "+".join(sorted(self.labels)) or BACKGROUND


def read_bed(path):
    for line in open(path):
        if line.strip() and not line.startswith(("#", "track", "browser")):
            c, s, e = line.split()[:3]
            yield c, int(s), int(e)


def build_atoms(beds, contigs, warn=print):
    """beds: {gene_set: bed_path}; contigs: [(name, length)] -> atoms tiling every contig.
    BED lines on contigs missing from `contigs` are ignored (e.g. "17" vs "chr17")."""
    ev, length = {c: [] for c, _ in contigs}, dict(contigs)
    for gs, path in beds.items():
        for c, s, e in read_bed(path):
            if c not in ev:
                warn(f"warning: {path}: contig {c} not in BAM, line ignored")
                continue
            s, e = max(0, s), min(e, length[c])
            if s < e:
                ev[c] += [(s, gs, 1), (e, gs, -1)]
    atoms = []
    for c, length in contigs:
        cuts = sorted({0, length, *(p for p, _, _ in ev[c])})
        evs, i, cnt = sorted(ev[c]), 0, {}
        for a, b in zip(cuts, cuts[1:]):
            while i < len(evs) and evs[i][0] == a:
                cnt[evs[i][1]] = cnt.get(evs[i][1], 0) + evs[i][2]
                i += 1
            lab = frozenset(g for g, n in cnt.items() if n > 0)
            if atoms and atoms[-1].contig == c and atoms[-1].labels == lab:
                atoms[-1] = Atom(c, atoms[-1].start, b, lab)
            else:
                atoms.append(Atom(c, a, b, lab))
    return atoms


def readable(gene_sets, classes, background=False):
    """Classes a reader entitled to `gene_sets` may open."""
    gs = set(gene_sets)
    return {c for c in classes if (c == BACKGROUND and background)
            or (c != BACKGROUND and ("*" in gs or gs & set(c.split("+"))))}
