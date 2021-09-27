#include "wright_fisher.hpp"

using namespace TCLAP;
using namespace std;

int main(int argc, char *argv[]) {
    CmdLine cmd{"neutral", ' ', "0.1"};
    OutputArgParse args(cmd);
    GenomeStructureArgParse args_genome(cmd);
    cmd.parse(argc, argv);

    generator.seed(args.seed.getValue());
    GenomeStructure genome = args_genome.get_model();
    NeutralModel neutral_fitness{};
    Population population(args.population_size.getValue(), genome, neutral_fitness);
    population.run(args.number_of_generations.getValue());
    population.trace.write_tsv(args.output_path.getValue());
    return 0;
}