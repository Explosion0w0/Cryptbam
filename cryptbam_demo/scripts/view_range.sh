#!/usr/bin/env bash
# Count the reads of a BAM in one demo gene region.
#   bash view_range.sh decrypted.bam BRCA1
#   bash view_range.sh decrypted.bam all     # every gene below
#   bash view_range.sh decrypted.bam         # list the genes
# Coordinates: GRCh38, 0-based BED as in generanges/*.bed (Ensembl canonical transcripts).
set -eu

declare -A REGION=(
  [BRCA1]="chr17:43044294-43125364"         # breast_cancer, ovarian_cancer
  [BRCA2]="chr13:32315507-32400268"         # breast_cancer, ovarian_cancer
  [PALB2]="chr16:23603164-23641310"         # breast_cancer
  [RAD51C]="chr17:58692601-58735611"        # ovarian_cancer
  [ALPHA-GLOBIN]="chr16:140000-200000"      # a-thalassemia (HBZ-HBA2-HBA1-HBQ1)
)
ORDER="BRCA1 BRCA2 PALB2 RAD51C ALPHA-GLOBIN"

if [ $# -lt 1 ] || [ ! -f "$1" ]; then
  echo "usage: bash $0 <file.bam> [gene|all]" >&2
  [ $# -ge 1 ] && echo "error: file not found: $1" >&2
  exit 1
fi
BAM=$1
command -v samtools >/dev/null || { echo "error: samtools not found (conda activate cryptbam)" >&2; exit 1; }
if [ $# -lt 2 ]; then
  echo "genes:"; for g in $ORDER; do printf '  %-13s %s\n' "$g" "${REGION[$g]}"; done; exit 0
fi

GENE=$(echo "$2" | tr '[:lower:]' '[:upper:]')
case $GENE in ALPHA|HBA|HBA1|HBA2|A-THALASSEMIA) GENE=ALPHA-GLOBIN ;; esac
if [ "$GENE" = ALL ]; then GENES=$ORDER
elif [ -n "${REGION[$GENE]+x}" ]; then GENES=$GENE
else echo "error: unknown gene '$2'; choose from: $ORDER all" >&2; exit 1
fi

# region queries need an index; build one next to the BAM if it is missing
[ -f "$BAM.bai" ] || [ -f "${BAM%.bam}.bai" ] || { echo "(indexing $BAM)" >&2; samtools index "$BAM"; }

# samtools uses 1-based inclusive regions; the table above is 0-based half-open
printf '%-13s %-26s %10s\n' gene region reads
for g in $GENES; do
  IFS=':-' read -r c s e <<< "${REGION[$g]}"
  n=$(samtools view -c "$BAM" "$c:$((s + 1))-$e")
  printf '%-13s %-26s %10s%s\n' "$g" "$c:$((s + 1))-$e" "$n" "$([ "$n" = 0 ] && echo '   (no reads: no access or not covered)')"
done
