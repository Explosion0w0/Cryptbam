#!/usr/bin/env bash
# Prepare the demo outputs for IGV (run after scripts/demo_real.sh):
#  * decrypted depth files -> .bedgraph without the "#" alignment padding lines
#  * demo/igv/*.batch      -> IGV Desktop: File > Run Batch Script (absolute paths for this machine)
#  * demo/igv/index.html   -> igv.js in the browser: npx http-server demo -p 8765, open /igv/
set -eu
cd "$(dirname "$0")/../demo"
for f in out/*.bdg; do
  b=$(basename "$f" .bdg)
  grep -v '^#' "$f" | sed "s/name=true_depth/name=\"$b true depth\"/" > "out/$b.bedgraph"
done
ROOT=$( (command -v wslpath >/dev/null && wslpath -m .) || pwd -W 2>/dev/null || pwd )
mkdir -p igv
batch() {  # file locus tracks...
  local out=igv/$1 locus=$2; shift 2
  { echo new; echo genome hg38
    for t in "$@"; do echo "load $ROOT/$t"; done
    echo "goto $locus"; } > "$out"
  echo "wrote $out"
}
batch alpha_thal.batch chr16:140000-200000 beds/a-thalassemia.bed   data/HG00684.alpha.bam out/HG00684.thal_analyst.bam out/HG00684.cnv_analyst.bedgraph   data/NA12878.genes.bam out/thal_analyst.bam out/cnv_analyst.bedgraph
batch access.batch chr17:43044294-43125364 beds/breast_cancer.bed beds/ovarian_cancer.bed   data/NA12878.genes.bam out/breast_analyst.bam out/ovarian_analyst.bam out/thal_analyst.bam out/cnv_analyst.bedgraph
