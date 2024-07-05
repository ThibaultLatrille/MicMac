CPU=8
for EXPERIMENT in ./config/*.yaml; do
python3 simulated_experiment.py -c $(basename "${EXPERIMENT}") -j ${CPU}
done
