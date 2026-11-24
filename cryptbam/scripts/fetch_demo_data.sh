#!/usr/bin/env bash
# Download a small real BAM for the demo: NA12878 (1000 Genomes 30x, GRCh38, public,
# no access request) restricted to the demo gene sets +-5 kb (breast_cancer,
# ovarian_cancer, a-thalassemia), plus HG00684, an alpha-thalassemia carrier. Only the needed CRAM
# slices and reference bases are fetched over HTTP (needs samtools with libcurl).
set -euo pipefail
cd "$(dirname "$0")/../demo"
CRAM=https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram
REF=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
mkdir -p data
FLANK=5000
regions=$(cat beds/*.bed | sort -k1,1 -k2,2n -u | awk -v f=$FLANK '{s=$2-f; if(s<0)s=0; print $1":"s+1"-"$3+f}')
echo "regions: $regions"
[ -f data/NA12878.genes.bam ] ||   # -M: multi-region iterator, each read once
  samtools view -M -T "$REF" -b -o data/NA12878.genes.bam "$CRAM" $regions
samtools index data/NA12878.genes.bam
# reference for the same windows (for bcftools / IGV), sequence names = chromosome names
if [ ! -f data/ref.fa.gz ]; then { python3 - "$REF" $regions <<'PY'
import sys, subprocess
ref, regs = sys.argv[1], sys.argv[2:]
for c in sorted({r.split(":")[0] for r in regs}, key=lambda c: int(c[3:])):
    L = int(next(l.split("\t")[1] for l in subprocess.run(["curl", "-s", ref + ".fai"], capture_output=True, text=True).stdout.splitlines() if l.split("\t")[0] == c))
    out, pos = [], 0
    for r in [r for r in regs if r.split(":")[0] == c]:
        s, e = map(int, r.split(":")[1].split("-"))
        fa = subprocess.run(["samtools", "faidx", ref, r], capture_output=True, text=True, check=True).stdout
        out.append("N" * (s - 1 - pos) + "".join(fa.split("\n")[1:]))
        pos = e
    out.append("N" * (L - pos))
    s = "".join(out)
    sys.stdout.write(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
PY
} | bgzip -c > data/ref.fa.gz.tmp
  mv data/ref.fa.gz.tmp data/ref.fa.gz && samtools faidx data/ref.fa.gz
fi
# a small plain-text SAM excerpt for the format-conversion demo: header trimmed to the
# primary chromosomes (the full GRCh38 header has 3366 @SQ lines) + 2000 BRCA1 reads
{ samtools view -H data/NA12878.genes.bam | grep -vP '^@SQ	SN:(?!chr([0-9]+|X|Y|M)	)'
  samtools view data/NA12878.genes.bam chr17:43044294-43125364 | head -n 2000 || true
} > data/NA12878.BRCA1.excerpt.sam

# Coverage-protection demo: HG00684 (Southern Han Chinese, 1000 Genomes 30x) carries a
# heterozygous --SEA alpha-thalassemia deletion (~20 kb, HBA1/HBA2); NA12878 above is normal there.
[ -f data/HG00684.alpha.bam ] || samtools view -M -T "$REF" -b -o data/HG00684.alpha.bam   https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241741/HG00684.final.cram chr16:130001-210000
samtools index data/HG00684.alpha.bam
rm -f ./*.crai                      # htslib caches remote CRAM indexes in the working directory
ls -la data
