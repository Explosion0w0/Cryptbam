#!/usr/bin/env bash
# Live demo on real public data (1000 Genomes 30x, GRCh38).
# Gene sets: breast_cancer (BRCA1, BRCA2, PALB2), ovarian_cancer (BRCA1, BRCA2, RAD51C),
#            a-thalassemia (alpha-globin locus, HBA1/HBA2).
#   bash scripts/fetch_demo_data.sh      # once, ~8 MB download (~10 min)
#   bash scripts/demo_real.sh            # PAUSE=1 bash scripts/demo_real.sh  -> wait for Enter between steps
set -eu          # no pipefail: demo commands pipe into head
cd "$(dirname "$0")/.."
export PYTHONPATH=$PWD
D=demo; O=$D/out; BAM=$D/data/NA12878.genes.bam; C4=$O/NA12878.bam.c4gh
BRCA1=chr17:43044294-43125364; ALPHA=chr16:140000-200000; ADEL=165000-185000
mkdir -p $O

step() { printf '\n\033[1;36m==== %s ====\033[0m\n' "$*"; [ "${PAUSE:-0}" = 1 ] && read -rp "(Enter) " _ || true; }
run()  { printf '\033[33m$ %s\033[0m\n' "$*"; eval "$@"; }

step "0. Keys: one key pair per identity (outsider = valid key, not in the policy)"
for n in breast_analyst ovarian_analyst thal_analyst cnv_analyst owner outsider; do
  [ -f $D/keys/$n.sec ] || python3 -m cryptbam keygen $n --outdir $D/keys
done
run "ls $D/keys"

step "1. A normal SAM file: SAM -> BAM -> SAM round trip with samtools"
S=$D/data/NA12878.BRCA1.excerpt.sam
run "grep -v '^@' $S | head -2 | cut -c1-150"
run "samtools view -b -o $O/excerpt.bam $S && samtools index $O/excerpt.bam"
run "samtools view -h -o $O/excerpt.roundtrip.sam $O/excerpt.bam"
run "diff <(grep -v '^@' $S) <(grep -v '^@' $O/excerpt.roundtrip.sam) && echo 'reads identical after SAM->BAM->SAM'"
run "samtools flagstat $O/excerpt.bam | head -2"

step "2. Plain BAM: whoever has the file sees everything"
run "samtools idxstats $BAM | awk '\$3>0'"
run "samtools view -H $BAM | grep -m1 '^@PG.*CL:' | cut -c1-160    # pipeline paths / sample IDs"
run "samtools view -c $BAM $ALPHA    # alpha-globin (thalassemia) reads, readable by anyone"

step "3. Standard Crypt4GH (official tool, one key for the whole file)"
run "crypt4gh encrypt --recipient_pk $D/keys/breast_analyst.pub < $BAM > $O/standard.bam.c4gh"
run "crypt4gh decrypt --sk $D/keys/outsider.sec < $O/standard.bam.c4gh > /dev/null || echo 'outsider: blocked'"
run "crypt4gh decrypt --sk $D/keys/breast_analyst.sec < $O/standard.bam.c4gh > $O/standard.dec.bam && samtools index $O/standard.dec.bam"
run "samtools view -c $O/standard.dec.bam $ALPHA    # breast analyst still reads alpha-globin: all or nothing"

step "4. cryptbam: one data key per access class, wrapped per reader"
run "cat $D/beds/*.bed"
run "python3 -m cryptbam encrypt --bam $BAM --policy $D/policy.json --out $O/NA12878 > $O/encrypt.log; grep depth= $O/encrypt.log"
run "ls -l $O/NA12878.*"

step "5. Correct key / wrong key / tampered file"
run "python3 -m cryptbam decrypt --sk $D/keys/breast_analyst.sec --in $C4 --out $O/breast_analyst.bam"
run "crypt4gh decrypt --sk $D/keys/breast_analyst.sec < $C4 | cmp - $O/breast_analyst.bam && echo 'official crypt4gh: identical bytes'"
run "samtools index $O/breast_analyst.bam && samtools view -c $O/breast_analyst.bam $BRCA1"
run "samtools view -c $O/breast_analyst.bam $ALPHA    # 0: alpha-globin is not in the breast_cancer gene set"
run "python3 -m cryptbam decrypt --sk $D/keys/outsider.sec --in $C4 --out /dev/null || true"
run "crypt4gh decrypt --sk $D/keys/outsider.sec < $C4 > /dev/null || echo 'official tool: outsider blocked too'"
cp $C4 $O/tampered.bam.c4gh
run "printf '\\xff' | dd of=$O/tampered.bam.c4gh bs=1 seek=\$((\$(stat -c%s $C4) - 5)) conv=notrunc 2>/dev/null    # flip 1 byte in the last segment"
run "python3 -m cryptbam decrypt --sk $D/keys/breast_analyst.sec --in $O/tampered.bam.c4gh --out /dev/null || true"

step "6. Same file, different identities"
for n in ovarian_analyst thal_analyst owner; do
  python3 -m cryptbam decrypt --sk $D/keys/$n.sec --in $C4 --out $O/$n.bam 2>/dev/null
  crypt4gh decrypt --sk $D/keys/$n.sec < $C4 | cmp -s - $O/$n.bam && samtools index $O/$n.bam
done
run "python3 scripts/identities.py $D"
run "python3 -m cryptbam inspect --sk $D/keys/thal_analyst.sec --in $C4    # audit: segments this key can open"

step "7. Coverage: reads have flat depth; true depth is a separate encrypted bedGraph"
run "samtools depth -a -r $BRCA1 $BAM | awk '{s+=\$3}END{printf \"plain BAM   BRCA1 mean depth %.1f\\n\", s/NR}'"
run "samtools depth -a -r $BRCA1 $O/breast_analyst.bam | awk '{s+=\$3}END{printf \"breast view BRCA1 mean depth %.1f\\n\", s/NR}'"
run "python3 -m cryptbam decrypt --sk $D/keys/cnv_analyst.sec --in $O/NA12878.depth.bdg.c4gh --out $O/cnv_analyst.bdg"
run "grep -v '^#' $O/cnv_analyst.bdg | cut -f1 | uniq -c    # only the three gene sets, no background"
run "python3 -m cryptbam decrypt --sk $D/keys/breast_analyst.sec --in $O/NA12878.depth.bdg.c4gh --out /dev/null || true"

step "8. Utility: variant calling in BRCA1, plain BAM vs breast analyst's view (K=10)"
call() { bcftools mpileup -f $D/data/ref.fa.gz -r $BRCA1 "$1" 2>/dev/null | bcftools call -mv 2>/dev/null |
         grep -v '^#' | awk '{split($10,g,":"); print $1,$2,$4,$5,g[1]}' | sort; }
call $BAM > $O/plain.calls; call $O/breast_analyst.bam > $O/breast.calls
run "echo \"plain: \$(wc -l < $O/plain.calls) calls, breast view: \$(wc -l < $O/breast.calls), shared: \$(comm -12 $O/plain.calls $O/breast.calls | wc -l)\""
run "comm -3 $O/plain.calls $O/breast.calls | cut -c1-80 | head -6    # differences = STR / homopolymer indels"

step "9. Coverage protection: alpha-thalassemia carrier (HG00684) vs normal (NA12878)"
# HG00684 carries a heterozygous --SEA deletion (~20 kb, HBA1/HBA2). Same policy, same keys.
python3 -m cryptbam encrypt --bam $D/data/HG00684.alpha.bam --policy $D/policy.json --out $O/HG00684 > /dev/null
python3 -m cryptbam decrypt --sk $D/keys/thal_analyst.sec --in $O/HG00684.bam.c4gh --out $O/HG00684.thal_analyst.bam 2>/dev/null
python3 -m cryptbam decrypt --sk $D/keys/cnv_analyst.sec --in $O/HG00684.depth.bdg.c4gh --out $O/HG00684.cnv_analyst.bdg 2>/dev/null
samtools index $O/HG00684.thal_analyst.bam
echo "thal_analyst may read the alpha-globin reads (e.g. point mutations) but has no coverage key:"
echo "can it tell who is a carrier from the reads alone?"
run "python3 scripts/depth_bins.py $ALPHA 5000 $ADEL \
  'HG00684 plain=$D/data/HG00684.alpha.bam' 'NA12878 plain=$BAM' \
  'HG00684 thal=$O/HG00684.thal_analyst.bam' 'NA12878 thal=$O/thal_analyst.bam'"
echo "cnv_analyst holds the coverage key: the true depth shows the deletion"
run "python3 scripts/depth_bins.py $ALPHA 5000 $ADEL 'HG00684 cnv.bdg=$O/HG00684.cnv_analyst.bdg' 'NA12878 cnv.bdg=$O/cnv_analyst.bdg'"
run "python3 -m cryptbam decrypt --sk $D/keys/thal_analyst.sec --in $O/HG00684.depth.bdg.c4gh --out /dev/null || true"
