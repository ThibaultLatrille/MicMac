CPU=8
for EXPERIMENT in ./config/constant_pop_size_many_sites*.yaml; do
python3 simulated_experiment.py -c $(basename "${EXPERIMENT}") -j ${CPU}
done
