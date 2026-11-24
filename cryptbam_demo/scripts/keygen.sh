mkdir -p ./keys
for n in breast_analyst ovarian_analyst thal_analyst cnv_analyst owner outsider; do
  [ -f ./keys/$n.sec ] || python3 -m cryptbam keygen $n --outdir ./keys
done
