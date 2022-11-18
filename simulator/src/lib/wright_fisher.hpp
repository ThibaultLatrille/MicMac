#include <chrono>
#include <map>
#include <set>
#include "argparse.hpp"
#include "fitness.hpp"
#include "random.hpp"
#include "statistic.hpp"
#include "tree.hpp"

struct TimeElapsed {
    double mutation{0.0};
    double selection{0.0};
    double fixation{0.0};
    double summary_stats{0.0};
};

TimeElapsed timer;

using namespace TCLAP;
using namespace std;

#define duration(a) chrono::duration_cast<chrono::nanoseconds>(a).count()
#define timeNow() chrono::high_resolution_clock::now()

typedef chrono::high_resolution_clock::time_point TimeVar;


class GenomeStructure : public vector<double> {
  public:
    u_long number_loci{};
    double mutation_mean_effect_size{};
    double mutation_rate_per_loci_per_generation{};
    double background_genotype{0};
    double expected_var{0};
    double std_environment{0};

    GenomeStructure() = default;
    ~GenomeStructure() = default;
    explicit GenomeStructure(u_long n, double radius, double mu, double variance_environment)
        : vector<double>(n, 0),
          number_loci{n},
          mutation_mean_effect_size{radius},
          mutation_rate_per_loci_per_generation{mu},
          std_environment{sqrt(variance_environment)} {
        for (u_long i = 0; i < number_loci; ++i) {
            (*this)[i] = mutation_mean_effect_size * normal_distrib(generator);
        }
        expected_var = 2 * mutation_rate_per_loci_per_generation * sum_squarred(*this);
    }

    double expected_mutations() const {
        return 2 * static_cast<double>(number_loci) * mutation_rate_per_loci_per_generation;
    };

    double expected_variance() const { return expected_var; };
};

class GenomeStructureArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    ~GenomeStructureArgParse() = default;
    explicit GenomeStructureArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}

    TCLAP::ValueArg<u_long> number_loci{
        "", "number_loci", "Number of loci", false, 10000, "u_long", cmd};
    TCLAP::ValueArg<double> mutation_mean_effect_size{
        "", "mutation_mean_effect_size", "Mutation mean effect size", false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> mutation_rate_per_loci_per_generation{"",
        "mutation_rate_per_loci_per_generation", "Mutation rate per locus per generation", false,
        1e-5, "double", cmd};
    TCLAP::ValueArg<double> variance_environment{
        "", "variance_environment", "Environmental variance (Ve).", false, 0.0, "double", cmd};

    GenomeStructure get_model() {
        return GenomeStructure(number_loci.getValue(), mutation_mean_effect_size.getValue(),
            mutation_rate_per_loci_per_generation.getValue(), variance_environment.getValue());
    }
};

class Individual {
  public:
    // If the value of the map is false, the locus is heterozygous.
    // If the value of the map is true, the locus is homozygous for the derived allele.
    unordered_map<u_long, bool> loci{};
    double genotype{0.0};
    double environment{0.0};
    double phenotype{0.0};
    double fitness{1.0};

    Individual() = default;
    ~Individual() = default;
    Individual(Individual const &mum, Individual const &dad, GenomeStructure const &genome,
        FitnessModel const &fitness_model)
        : loci{} {
        parent_loci_sampling(mum);
        parent_loci_sampling(dad);
        genotype = compute_genotype(genome);
        if (genome.std_environment != 0.0) {
            environment = genome.std_environment * normal_distrib(generator);
        }
        phenotype = genotype + environment;
        update_fitness(fitness_model);
    };

    double heterozygosity(GenomeStructure const &genome) const {
        u_long h{0};
        for (auto const &locus : loci) {
            if (!locus.second) { h++; }
        }
        return static_cast<double>(h) / static_cast<double>(genome.number_loci);
    }

    double compute_genotype(GenomeStructure const &genome) const {
        double p = genome.background_genotype;
        for (auto const &locus : loci) { p += genome[locus.first] * (1 + locus.second); }
        return p;
    }

    void parent_loci_sampling(Individual const &parent) {
        for (auto const &locus_pair : parent.loci) {
            if (locus_pair.second or bernouilli_distrib(generator)) {
                auto it = loci.find(locus_pair.first);
                if (it != loci.end()) {
                    it->second = true;
                } else {
                    loci[locus_pair.first] = false;
                }
            }
        }
    }

    void update_fitness(FitnessModel const &fitness_model) {
        fitness = fitness_model.fitness(phenotype);
    };

    void mutation(GenomeStructure const &genome, FitnessModel const &fitness_model) {
        poisson_distribution<u_long> poisson_distrib(genome.expected_mutations());
        u_long poisson_draw = poisson_distrib(generator);

        if (poisson_draw > 0) {
            uniform_int_distribution<u_long> u_distrib(0, genome.number_loci - 1);

            for (u_long mutation{0}; mutation < poisson_draw; mutation++) {
                auto locus = u_distrib(generator);
                auto it = loci.find(locus);
                if (it != loci.end()) {
                    if (it->second) {
                        it->second = false;
                    } else {
                        loci.erase(it);
                    }
                    genotype -= genome.at(locus);
                } else {
                    loci[locus] = false;
                    genotype += genome.at(locus);
                }
                phenotype = genotype + environment;
            }
        }
        update_fitness(fitness_model);
    };
};


class PopulationSizeProcess {
  private:
    u_long population_size;
    u_long population_size_min;

    double brownian;
    double brownian_min;
    double brownian_sigma;

    double ou{0};
    double ou_sigma;
    double ou_theta;
    normal_distribution<double> normal_distrib;

  public:
    explicit PopulationSizeProcess(u_long population_size, u_long population_size_min,
        double brownian_sigma, double ou_sigma, double ou_theta)
        : population_size{population_size},
          population_size_min{population_size_min},
          brownian_sigma{brownian_sigma},
          ou_sigma{ou_sigma},
          ou_theta{ou_theta} {
        brownian = log(population_size);
        brownian_min = log(population_size_min);
        normal_distrib = normal_distribution<double>(0.0, 1.0);
    }

    void update() {
        ou += ou_sigma * normal_distrib(generator_pop_size) - ou_theta * ou;
        brownian += brownian_sigma * normal_distrib(generator_pop_size);
        if (brownian < brownian_min) { brownian = 2 * brownian_min - brownian; }
        population_size = static_cast<u_long>(exp(ou + brownian));
    }
    u_long get_population_size() const { return max(population_size, population_size_min); }
};


class PopulationSizeArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    ~PopulationSizeArgParse() = default;
    explicit PopulationSizeArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}

    TCLAP::ValueArg<u_long> population_size{
        "", "population_size", "Number of individuals", false, 100, "u_long", cmd};
    TCLAP::ValueArg<u_long> population_size_min{
        "", "population_size_min", "Minimum population size", false, 20, "u_long", cmd};
    TCLAP::ValueArg<double> brownian_sigma{"", "brownian_sigma",
        "The Brownian sigma (0<sigma) applied to Ne at each generation", false, 0.05, "double",
        cmd};
    TCLAP::ValueArg<double> noise_sigma{"", "noise_sigma",
        "The Ornstein–Uhlenbeck sigma (0<sigma) applied to Ne at each generation", false, 0.1,
        "double", cmd};
    TCLAP::ValueArg<double> noise_theta{"", "noise_theta",
        "The Ornstein–Uhlenbeck theta (0<=theta<1) applied to Ne at each generation", false, 0.9,
        "double", cmd};

    PopulationSizeProcess get_model() {
        return PopulationSizeProcess(population_size.getValue(), population_size_min.getValue(),
            brownian_sigma.getValue(), noise_sigma.getValue(), noise_theta.getValue());
    }
};

class Population {
  public:
    u_long elapsed{0};
    string name{"ROOT"};
    GenomeStructure genome;
    FitnessModel &fitness_function;
    PopulationSizeProcess pop_size;
    FitnessState fitness_state;
    mutable double heritability{0};

    // Individuals composing the population
    vector<Individual> parents{};
    unsigned nbr_fixations{0};

    ~Population() = default;
    explicit Population(
        const GenomeStructure &genome, FitnessModel &f_model, const PopulationSizeProcess &pop_size)
        : genome(genome), fitness_function{f_model}, pop_size(pop_size) {
        fitness_state = fitness_function.get_state();
        parents.resize(pop_size.get_population_size());
        selection_and_random_mating();
        mutation();
    };

    void update(Population const &pop) {
        elapsed = pop.elapsed;
        genome = pop.genome;
        pop_size = pop.pop_size;
        parents = pop.parents;
        fitness_state = pop.fitness_state;
        fitness_function.set_state(fitness_state);
        nbr_fixations = 0;
    }

    bool check() const {
        return all_of(parents.begin(), parents.end(), [this](Individual const &p) {
            return abs(p.genotype - p.compute_genotype(genome)) < 1e-3;
        });
    }

    void mutation() {
        for (auto &p : parents) { p.mutation(genome, fitness_function); }
    };

    vector<Individual> sampling_individuals(bool fitness) const {
        vector<double> fitness_vector(parents.size(), 0);
        transform(parents.begin(), parents.end(), fitness_vector.begin(),
            [&fitness](Individual const &b) { return (fitness ? b.fitness : 1.0); });
        discrete_distribution<u_long> parent_distrib(fitness_vector.begin(), fitness_vector.end());

        vector<double> x, y;
        vector<Individual> children;
        children.reserve(pop_size.get_population_size());
        for (u_long i{0}; i < pop_size.get_population_size(); i++) {
            u_long mum = parent_distrib(generator);
            u_long dad = parent_distrib(generator);
            x.push_back((parents[mum].phenotype + parents[dad].phenotype) / 2);
            children.emplace_back(parents[mum], parents[dad], genome, fitness_function);
            y.push_back(children.back().phenotype);
        }
        heritability = slope(x, y);
        return children;
    }

    void selection_and_random_mating() { parents = sampling_individuals(true); };

    unordered_map<u_long, u_long> site_frequency_spectrum() const {
        unordered_map<u_long, u_long> frequencies{};
        for (Individual const &p : parents) {
            for (pair<const u_long, bool> const &locus : p.loci) {
                frequencies[locus.first] += 1 + locus.second;
            }
        }
        unordered_map<u_long, u_long> sfs{};
        for (pair<const u_long, u_long> const &freq : frequencies) { sfs[freq.second] += 1; }
        return sfs;
    }

    tuple<double, double, double> theta(unordered_map<u_long, u_long> const &sfs) const {
        Stat watterson{}, tajima{}, fay_wu{};
        u_long n = pop_size.get_population_size() * 2;
        assert(parents.size() * 2 == n);
        for (u_long i{1}; i < n; i++) {
            double epsilon = {0.0};
            auto it = sfs.find(i);
            if (it != sfs.end()) {
                epsilon = static_cast<double>(it->second);
                assert(epsilon > 0);
            }
            double theta = static_cast<double>(i) * epsilon;
            watterson.add(theta, 1.0 / static_cast<double>(i));
            tajima.add(theta, static_cast<double>(n - i));
            fay_wu.add(theta, static_cast<double>(i));
        }
        double theta_watterson = watterson.mean() / static_cast<double>(genome.number_loci);
        double theta_tajima = tajima.mean() / static_cast<double>(genome.number_loci);
        double theta_fay_wu = fay_wu.mean() / static_cast<double>(genome.number_loci);
        return {theta_watterson, theta_tajima, theta_fay_wu};
    }

    double sampled_theta_pairwise() const {
        double h{0.0};
        auto kids = sampling_individuals(false);
        for (auto const &individual : kids) { h += individual.heterozygosity(genome); }
        return h / static_cast<double>(kids.size());
    }

    void fixation() {
        set<u_long> loci_tmp{};
        for (auto const &locus_pair : parents[0].loci) {
            if (locus_pair.second) { loci_tmp.emplace(locus_pair.first); }
        }
        for (u_long i{1}; i < parents.size(); i++) {
            set<u_long> loci_to_remove{};
            for (auto const &locus : loci_tmp) {
                auto it = parents[i].loci.find(locus);
                if (it == parents[i].loci.end() or !it->second) { loci_to_remove.emplace(locus); }
            }
            for (auto const &v : loci_to_remove) { loci_tmp.erase(v); }
            if (loci_tmp.empty()) { break; }
        }
        for (auto const &fixation : loci_tmp) {
            nbr_fixations++;
            genome.background_genotype += genome[fixation] * 2;
            genome[fixation] = -genome[fixation];
            for (auto &p : parents) { p.loci.erase(fixation); }
        }
        assert(check());
    }

    unordered_map<string, double> summary_states() const {
        unordered_map<string, double> stats{};

        vector<double> traits(parents.size(), 0);
        transform(parents.begin(), parents.end(), traits.begin(),
            [](Individual const &b) { return b.phenotype; });
        stats["Phenotype mean"] = mean(traits);
        stats["Phenotype var"] = variance(traits);

        transform(parents.begin(), parents.end(), traits.begin(),
            [](Individual const &b) { return b.genotype; });
        stats["Genotype mean"] = mean(traits);
        stats["Genotype var"] = variance(traits);

        transform(parents.begin(), parents.end(), traits.begin(),
            [](Individual const &b) { return b.environment; });
        stats["Environment var"] = variance(traits);

        transform(parents.begin(), parents.end(), traits.begin(),
            [](Individual const &b) { return b.fitness; });
        stats["Fitness mean"] = mean(traits);
        stats["Fitness var"] = variance(traits);
        return stats;
    }

    void save_node(Tree &tree, Tree::NodeIndex node) {
        tree.set_tag(node, "Population size", to_string(pop_size.get_population_size()));
        tree.set_tag(node, "Number of loci", to_string(genome.number_loci));
        tree.set_tag(
            node, "Mutational rate", to_string(genome.mutation_rate_per_loci_per_generation));
        tree.set_tag(node, "Mutational var", to_string(genome.expected_variance()));
        tree.set_tag(node, "Heritability", to_string(heritability));
        auto stats = summary_states();
        for (auto const &stat : stats) { tree.set_tag(node, stat.first, to_string(stat.second)); }

        fixation();
        auto sfs = site_frequency_spectrum();
        auto theta_tuple = theta(sfs);
        tree.set_tag(node, "Theta Watterson", to_string(get<0>(theta_tuple)));
        tree.set_tag(node, "Theta Tajima", to_string(get<1>(theta_tuple)));
        tree.set_tag(node, "Theta Fay Wu", to_string(get<2>(theta_tuple)));
        tree.set_tag(node, "Theta", to_string(sampled_theta_pairwise()));

        double d = nbr_fixations / static_cast<double>(genome.number_loci);
        tree.set_tag(node, "Fixations", to_string(nbr_fixations));
        tree.set_tag(node, "d", to_string(d));
        tree.set_tag(node, "q", to_string(d / tree.node_length(node)));
    }

    void run(u_long nbr_generations, Trace &trace) {
        cout << "Population " << name << " for " << nbr_generations << " generations." << endl;
        // Run under selection
        u_long last_pct = 0;
        for (u_long sample{1}; sample <= nbr_generations; sample++) {
            elapsed++;
            pop_size.update();
            trace.add("Lineage", name);
            trace.add("Generation", elapsed);
            trace.add("Population size", pop_size.get_population_size());

            TimeVar t_start = timeNow();
            selection_and_random_mating();
            trace.add("Heritability", heritability);
            timer.selection += duration(timeNow() - t_start);

            t_start = timeNow();
            trace.add("Mutational expected var", genome.expected_variance());
            vector<double> traits(parents.size(), 0);
            transform(parents.begin(), parents.end(), traits.begin(),
                [](Individual const &b) { return b.phenotype; });
            double var_mutational = -variance(traits);
            mutation();
            transform(parents.begin(), parents.end(), traits.begin(),
                [](Individual const &b) { return b.phenotype; });
            var_mutational += variance(traits);
            trace.add("Mutational var", var_mutational);
            timer.mutation += duration(timeNow() - t_start);

            t_start = timeNow();
            auto stats = summary_states();
            for (auto const &stat : stats) { trace.add(stat.first, stat.second); }
            timer.summary_stats += duration(timeNow() - t_start);

            u_long pct = 25 * sample / nbr_generations;
            fitness_function.update();
            if (pct != last_pct) {
                last_pct = pct;
                t_start = timeNow();
                fixation();
                timer.fixation += duration(timeNow() - t_start);
            }
        }
        fitness_state = fitness_function.get_state();
    };
};

class Process {
  private:
    static double generations_computed;
    Tree &tree;
    vector<Population> populations;

  public:
    explicit Process(Tree &tree, Population const &pop, Trace &trace) : tree{tree} {
        for (Tree::NodeIndex i = 0; i < static_cast<int>(tree.nb_nodes()); ++i) {
            populations.emplace_back(pop);
            populations.back().name = tree.node_name(i);
        }
        run_recursive(tree.root(), trace);
        timer_count();
    }


    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, Trace &trace) {
        // Substitutions of the DNA sequence is generated.

        if (!tree.is_root(node)) {
            populations.at(node).run(static_cast<u_long>(tree.node_length(node)), trace);
            generations_computed += tree.node_length(node);
            cout << generations_computed << " generations computed in total ("
                 << static_cast<int>(100 * generations_computed / tree.total_length())
                 << "%) at node " << tree.node_name(node) << " ("
                 << static_cast<int>(100 * tree.node_length(node) / tree.total_length()) << "%)."
                 << endl;
            populations.at(node).save_node(tree, node);
        }

        // Iterate through the direct children.
        for (auto &child : tree.children(node)) {
            populations.at(child).update(populations.at(node));
            run_recursive(child, trace);
        }
    }

    static void timer_count() {
        double total_time = timer.mutation + timer.selection + timer.fixation + timer.summary_stats;
        cout << setprecision(3) << total_time / 1e9 << "s total time" << endl;
        cout << 100 * timer.mutation / total_time << "% of time spent in mutation ("
             << timer.mutation / 1e9 << "s)" << endl;
        cout << 100 * timer.selection / total_time << "% of time spent in selection ("
             << timer.selection / 1e9 << "s)" << endl;
        cout << 100 * timer.fixation / total_time << "% of time spent in fixation ("
             << timer.fixation / 1e9 << "s)" << endl;
        cout << 100 * timer.summary_stats / total_time << "% of time spent in summary stats ("
             << timer.summary_stats / 1e9 << "s)" << endl;
    }
};

// Initialize static variables
double Process::generations_computed = 0.0;