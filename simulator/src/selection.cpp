#include "wright_fisher.hpp"

using namespace TCLAP;

int main(int argc, char *argv[]) {
    CmdLine cmd{"stabilizing", ' ', "0.1"};
    OutputArgParse args(cmd);
    GenomeStructureArgParse args_genome(cmd);
    PopulationSizeArgParse args_pop_size(cmd);
    OrnsteinUhlenbeckBiasArgParse args_fitness(cmd);
    MutationRateArgParse args_mutation_rate(cmd);
    cmd.parse(argc, argv);

    generator.seed(args.seed.getValue());
    generator_pop_size.seed(args.seed_pop_size.getValue());
    generator_mut_rate.seed(args.seed_mut_rate.getValue());
    OrnsteinUhlenbeckBiasModel fitness = args_fitness.get_model();
    Population population(args_genome.get_model(), fitness, args_pop_size.get_model(),
        args_mutation_rate.get_model());

    double root_age{args.number_of_generations.getValue()};
    Tree tree = args.tree_path.getValue().empty()
                    ? Tree(args.number_of_lineages.getValue(), root_age)
                    : Tree(args.tree_path.getValue());
    tree.set_root_age(root_age);

    Trace trace(args.ouput_tsv.getValue());
    population.run(args.burn_in.getValue(), trace);
    Process process(tree, population, trace);
    trace.write_tsv(args.output_path.getValue());
    tree.write_newick(args.output_path.getValue());
    return 0;
}