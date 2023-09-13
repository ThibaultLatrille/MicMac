FOLDER=$(pwd)
for EXPERIMENT in $(pwd)/experiments/fluc*ing; do
  echo "Updating snakemake in: ${EXPERIMENT}"
  rm -rf "${EXPERIMENT}"/Snakefile
  ln $(pwd)/Snakefile "${EXPERIMENT}"/Snakefile
  cd "${EXPERIMENT}"
  snakemake -j 8 --touch
  cd "${FOLDER}"
done
