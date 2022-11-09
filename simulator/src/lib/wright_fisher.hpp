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
    map<u_long, bool> loci{};
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
    double heritability{0};

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

    void selection_and_random_mating() {
        vector<double> fitness_vector(parents.size(), 0);
        transform(parents.begin(), parents.end(), fitness_vector.begin(),
            [](Individual const &b) { return b.fitness; });
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
        parents = std::move(children);
    };

    unordered_map<u_long, u_long> site_frequency_spectrum() const {
        unordered_map<u_long, u_long> frequencies{};
        for (Individual const &p : parents) {
            for (pair<const u_long, bool> const &locus : p.loci) {
                // check if the locus is in frequencies
                if (frequencies.find(locus.first) == frequencies.end()) {
                    frequencies[locus.first] = 1 + locus.second;
                } else {
                    frequencies[locus.first] += 1 + locus.second;
                };
            }
        }
        unordered_map<u_long, u_long> sfs{};
        for (pair<const u_long, u_long> const &freq : frequencies) {
            // check if the locus is in frequencies
            if (frequencies.find(freq.second) == frequencies.end()) {
                sfs[freq.second] = 1;
            } else {
                sfs[freq.second] += 1;
            }
        }
        return sfs;
    }

    double theta(unordered_map<u_long, u_long> const &sfs, string const &method) const {
        double tot_weight{0};
        double tot_theta{0};
        u_long n = pop_size.get_population_size() * 2;
        assert(parents.size() * 2 == n);
        for (u_long i{1}; i < n; i++) {
            double epsilon = {0.0};
            auto it = sfs.find(i);
            if (it != sfs.end()) {
                epsilon = static_cast<double>(it->second);
                assert(epsilon > 0);
            }
            assert(epsilon <= genome.number_loci);

            double weight{0.0};
            if (method == "watterson") {
                weight = 1.0 / static_cast<double>(i);
            } else if (method == "tajima") {
                weight = static_cast<double>(pop_size.get_population_size() * 2 - i);
            } else if (method == "fay_wu") {
                weight = static_cast<double>(i);
            } else {
                cerr << "Unknown method for theta: " << method << endl;
                exit(1);
            }
            tot_weight += weight;
            tot_theta += weight * static_cast<double>(i) * epsilon;
        }
        return tot_theta / (tot_weight * static_cast<double>(genome.number_loci));
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
            timer.selection += duration(timeNow() - t_start);
            trace.add("Heritability", heritability);

            trace.add("Mutational expected var", genome.expected_variance());
            vector<double> traits(parents.size(), 0);
            transform(parents.begin(), parents.end(), traits.begin(),
                [](Individual const &b) { return b.phenotype; });

            double var_mutational = -variance(traits);
            t_start = timeNow();
            mutation();
            timer.mutation += duration(timeNow() - t_start);

            t_start = timeNow();
            transform(parents.begin(), parents.end(), traits.begin(),
                [](Individual const &b) { return b.phenotype; });
            var_mutational += variance(traits);
            trace.add("Mutational var", var_mutational);
            trace.add("Phenotype mean", mean(traits));
            trace.add("Phenotype var", variance(traits));

            transform(parents.begin(), parents.end(), traits.begin(),
                [](Individual const &b) { return b.genotype; });
            trace.add("Genotype mean", mean(traits));
            trace.add("Genotype var", variance(traits));

            transform(parents.begin(), parents.end(), traits.begin(),
                [](Individual const &b) { return b.environment; });
            trace.add("Environment var", variance(traits));

            transform(parents.begin(), parents.end(), traits.begin(),
                [](Individual const &b) { return b.fitness; });
            trace.add("Fitness mean", mean(traits));
            trace.add("Fitness var", variance(traits));

            auto sfs = site_frequency_spectrum();
            trace.add("Theta Watterson", theta(sfs, "watterson"));
            trace.add("Theta Tajima", theta(sfs, "tajima"));
            trace.add("Theta Fay Wu", theta(sfs, "fay_wu"));
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
        cout << 100 * timer.summary_stats / total_time << "% of time spent in fixation ("
             << timer.summary_stats / 1e9 << "s)" << endl;
    }
};

// Initialize static variables
double Process::generations_computed = 0.0;