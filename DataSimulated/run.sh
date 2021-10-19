CPU=8
for EXPERIMENT in ./*.yaml; do
# shellcheck disable=SC2046
python3 ./simulated_experiment.py -c $(basename "${EXPERIMENT}") -j ${CPU}
done
