#include "argparse.hpp"
#include "fitness.hpp"
#include "wright_fisher.hpp"

using namespace TCLAP;
using namespace std;

int main(int argc, char *argv[]) {
    CmdLine cmd{"stabilizing", ' ', "0.1"};
    OutputArgParse args(cmd);
    GenomeStructureArgParse args_genome(cmd);
    MovingOptimumArgParse args_moving_optimum(cmd);
    cmd.parse(argc, argv);

    generator.seed(args.seed.getValue());
    GenomeStructure genome = args_genome.get_model();
    MovingOptimumModel moving_optimum_fitness = args_moving_optimum.get_moving_optimum_model(args_genome.number_loci.getValue());
    Population population(args.population_size.getValue(), genome, moving_optimum_fitness);

    Trace trace{};
    population.run(args.burn_in.getValue(), trace);
    Process process(population, args.number_of_generations.getValue(), args.number_of_lineages.getValue(), trace);
    trace.write_tsv(args.output_path.getValue());
    return 0;
}