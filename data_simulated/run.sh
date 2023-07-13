CPU=8
for EXPERIMENT in ./constant_pop_size*.yaml; do
python3 simulated_experiment.py -c $(basename "${EXPERIMENT}") -j ${CPU}
done
