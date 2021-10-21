#include "wright_fisher.hpp"

using namespace TCLAP;

int main(int argc, char *argv[]) {
    CmdLine cmd{"stabilizing", ' ', "0.1"};
    OutputArgParse args(cmd);
    GenomeStructureArgParse args_genome(cmd);
    PopulationSizeArgParse args_pop_size(cmd);
    GaussianArgParse args_gaussian(cmd);
    cmd.parse(argc, argv);

    generator.seed(args.seed.getValue());
    generator_pop_size.seed(args.seed_pop_size.getValue());
    GaussianModel gaussian_fitness = args_gaussian.get_gaussian_model();
    Population population(args_genome.get_model(), gaussian_fitness, args_pop_size.get_model());

    double root_age{args.number_of_generations.getValue()};
    Tree tree = args.newick_path.getValue().empty()
                    ? Tree(args.number_of_lineages.getValue(), root_age)
                    : Tree(args.newick_path.getValue());
    tree.set_root_age(root_age);

    Trace trace{};
    population.run(args.burn_in.getValue(), trace);
    Process process(tree, population, trace);
    trace.write_tsv(args.output_path.getValue());
    tree.write_newick(args.output_path.getValue());
    return 0;
}