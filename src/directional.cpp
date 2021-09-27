#include "argparse.hpp"
#include "fitness.hpp"
#include "wright_fisher.hpp"

using namespace TCLAP;
using namespace std;

int main(int argc, char *argv[]) {
    CmdLine cmd{"stabilizing", ' ', "0.1"};
    OutputArgParse args(cmd);
    GenomeStructureArgParse args_genome(cmd);
    cmd.parse(argc, argv);

    generator.seed(args.seed.getValue());
    GenomeStructure genome = args_genome.get_model();
    DirectionalModel directional_fitness{};
    Population population(args.population_size.getValue(), genome, directional_fitness);

    Trace trace{};
    population.run(args.burn_in.getValue(), trace);
    Process process(population, args.number_of_generations.getValue(), trace);
    trace.write_tsv(args.output_path.getValue());
    return 0;
}