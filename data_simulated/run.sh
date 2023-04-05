CPU=6
for EXPERIMENT in ./constant_pop_size_L*.yaml; do
python3 simulated_experiment.py -c $(basename "${EXPERIMENT}") -j ${CPU}
done
