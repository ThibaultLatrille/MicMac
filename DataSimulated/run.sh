CPU=8
for EXPERIMENT in ./*.yaml; do
python3 ./simulated_experiment.py -c $(basename "${EXPERIMENT}") -j ${CPU}
done
