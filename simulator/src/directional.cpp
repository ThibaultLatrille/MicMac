#include "wright_fisher.hpp"

using namespace TCLAP;
using namespace std;

int main(int argc, char *argv[]) {
    CmdLine cmd{"stabilizing", ' ', "0.1"};
    OutputArgParse args(cmd);
    GenomeStructureArgParse args_genome(cmd);
    PopulationSizeArgParse args_pop_size(cmd);
    DirectionalArgParse args_directional(cmd);
    cmd.parse(argc, argv);

    generator.seed(args.seed.getValue());
    DirectionalModel directional_fitness = args_directional.get_directional_model();
    Population population(args_genome.get_model(), directional_fitness, args_pop_size.get_model());

    Trace trace{};
    population.run(args.burn_in.getValue(), trace);
    Process process(population, args.number_of_generations.getValue(), args.number_of_lineages.getValue(), trace);
    trace.write_tsv(args.output_path.getValue());
    return 0;
}