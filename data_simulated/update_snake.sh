FOLDER=$(pwd)
for EXPERIMENT in $(pwd)/experiments/*; do
  echo "Updating snakemake in: ${EXPERIMENT}"
  rm -rf "${EXPERIMENT}"/Snakefile
  ln $(pwd)/Snakefile "${EXPERIMENT}"/Snakefile
  cd "${EXPERIMENT}"
  snakemake --touch -j 8 --rerun-incomplete
  cd "${FOLDER}"
done
