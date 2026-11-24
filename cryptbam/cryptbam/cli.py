"""python -m cryptbam keygen|encrypt|decrypt|inspect  (see --help of each)."""
import argparse
import json
import os
import sys

from . import c4gh
from .pipeline import encrypt, load_policy


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cryptbam")
    sub = ap.add_subparsers(dest="cmd", required=True)
    k = sub.add_parser("keygen", help="create NAME.sec / NAME.pub")
    k.add_argument("name")
    k.add_argument("--outdir", default="keys")
    e = sub.add_parser("encrypt", help="BAM -> PREFIX.bam.c4gh + PREFIX.depth.bdg.c4gh")
    e.add_argument("--bam", required=True)
    e.add_argument("--policy", required=True)
    e.add_argument("--out", required=True, help="output prefix")
    e.add_argument("--depth", type=int, default=10,
                   help="default scaling for gene sets without 'depth' in the policy (0 = real reads)")
    e.add_argument("--set-depth", action="append", default=[], metavar="GENESET=K",
                   help="override one gene set, e.g. --set-depth a-thalassemia=0 --set-depth breast_cancer=30")
    e.add_argument("--no-background", action="store_true", help="drop reads outside every BED")
    for name, hlp in (("decrypt", "decrypt what this key may read"),
                      ("inspect", "try every segment with this key (audit)")):
        d = sub.add_parser(name, help=hlp)
        d.add_argument("--sk", required=True)
        d.add_argument("--in", dest="inp", required=True)
        if name == "decrypt":
            d.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    if a.cmd == "keygen":
        os.makedirs(a.outdir, exist_ok=True)
        c4gh.keygen(f"{a.outdir}/{a.name}.sec", f"{a.outdir}/{a.name}.pub")
    elif a.cmd == "encrypt":
        over = {s.split("=")[0]: int(s.split("=")[1]) for s in a.set_depth}
        access, _ = encrypt(a.bam, load_policy(a.policy, over), a.out, a.depth,
                            background=not a.no_background)
        print(json.dumps(access, indent=1))
    else:
        with open(a.inp, "rb") as f, open(a.out if a.cmd == "decrypt" else os.devnull, "wb") as o:
            print(json.dumps(c4gh.decrypt(f, o, c4gh.load_sk(a.sk), a.cmd == "decrypt")), file=sys.stderr)


def run():
    try:
        main()
    except (PermissionError, ValueError, FileNotFoundError) as err:
        sys.exit(f"cryptbam: error: {err}")
