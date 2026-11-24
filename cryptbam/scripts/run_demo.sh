#!/usr/bin/env bash
# simulate -> keys -> encrypt -> decrypt per analyst -> samtools
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/simulate.py --outdir data
for n in breast_analyst ovarian_analyst thal_analyst cnv_analyst owner; do
  [ -f data/keys/$n.sec ] || python3 -m cryptbam keygen $n --outdir data/keys
done
# depth per gene set comes from data/policy.json; override here if needed
python3 -m cryptbam encrypt --bam data/sample.bam --policy data/policy.json --out out/demo/sample \
        --set-depth breast_cancer=20
for n in breast_analyst ovarian_analyst thal_analyst owner; do
  python3 -m cryptbam decrypt --sk data/keys/$n.sec --in out/demo/sample.bam.c4gh --out out/demo/$n.bam
  if command -v crypt4gh >/dev/null; then      # official GA4GH tool gives identical bytes
    crypt4gh decrypt --sk data/keys/$n.sec < out/demo/sample.bam.c4gh | cmp - out/demo/$n.bam && echo "$n: identical"
  fi
  samtools index out/demo/$n.bam
  echo "$n: $(samtools view -c out/demo/$n.bam) reads"
done
python3 -m cryptbam decrypt --sk data/keys/cnv_analyst.sec --in out/demo/sample.depth.bdg.c4gh \
        --out out/demo/cnv_analyst.depth.bdg
python3 scripts/compare.py      # plain BAM vs standard Crypt4GH vs cryptbam
