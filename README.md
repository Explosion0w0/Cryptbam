# cryptbam

Region-based encryption of NGS BAM files, compatible with GA4GH Crypt4GH. Each gene set gets its own data key, so one encrypted file gives each analyst only the regions they are allowed to read, and depth scaling hides copy-number changes.

| Folder | Contents |
|---|---|
| [`cryptbam/`](cryptbam/) | The tool (`python -m cryptbam keygen / encrypt / decrypt / inspect`), simulation and evaluation scripts, the real-data demo. See [cryptbam/README.md](cryptbam/README.md) and [cryptbam/demo/DEMO.md](cryptbam/demo/DEMO.md). |
| [`cryptbam_demo/`](cryptbam_demo/) | Minimal demo: one patient BAM (1000 Genomes HG00684, public), reference, three risk-gene BEDs, and [`demo.txt`](cryptbam_demo/demo.txt) to run step by step. |

Requires Linux or WSL (pysam does not support Windows): `pip install pysam cryptography numpy` plus `samtools`.

Private keys (`*.sec`, `keys/`) and generated outputs are not in the repo; `cryptbam_demo/scripts/keygen.sh` creates the demo keys.
