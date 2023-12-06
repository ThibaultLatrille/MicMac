# shellcheck disable=SC2164
# shellcheck disable=SC2046
FOLDER=$(pwd)
for EXPERIMENT in $(pwd)/experiments/*_u; do
  echo "Updating snakemake in: ${EXPERIMENT}"
  rm -rf "${EXPERIMENT}"/Snakefile
  ln $(pwd)/Snakefile "${EXPERIMENT}"/Snakefile
  cd "${EXPERIMENT}"
  snakemake -j 8 -k
  cd "${FOLDER}"
done
