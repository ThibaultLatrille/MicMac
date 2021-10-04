#include <chrono>
#include <map>
#include <set>
#include "argparse.hpp"
#include "fitness.hpp"
#include "random.hpp"
#include "statistic.hpp"

#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

typedef std::chrono::high_resolution_clock::time_point TimeVar;

using namespace TCLAP;
using namespace std;

struct TimeElapsed {
    double mutation{0.0};
    double selection{0.0};
    double fixation{0.0};
};

TimeElapsed timer;

class GenomeStructure : public vector<double> {
  public:
    u_long number_loci{};
    double mutation_mean_effect_size{};
    double mutation_rate_per_loci_per_generation{};
    double background_phenotype{0};

    GenomeStructure() = default;
    ~GenomeStructure() = default;
    explicit GenomeStructure(u_long n, double radius, double mu)
        : vector<double>(n, 0),
          number_loci{n},
          mutation_mean_effect_size{radius},
          mutation_rate_per_loci_per_generation{mu} {
        for (u_long i = 0; i < number_loci; ++i) {
            (*this)[i] = mutation_mean_effect_size * normal_distrib(generator);
        }
    }

    double expected_mutations() const {
        return 2 * static_cast<double>(number_loci) * mutation_rate_per_loci_per_generation;
    };

    double expected_variance() const {
        return expected_mutations() * mutation_mean_effect_size * mutation_mean_effect_size;
    };
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

    GenomeStructure get_model() {
        return GenomeStructure(number_loci.getValue(), mutation_mean_effect_size.getValue(),
            mutation_rate_per_loci_per_generation.getValue());
    }
};

class Individual {
  public:
    map<u_long, bool> loci{};
    double phenotype{0.0};
    double fitness{1.0};

    Individual() = default;
    ~Individual() = default;
    Individual(Individual const &mum, Individual const &dad, GenomeStructure const &genome,
        FitnessModel const &fitness_model)
        : loci{} {
        parent_loci_sampling(mum);
        parent_loci_sampling(dad);
        phenotype = compute_phenotype(genome);
        update_fitness(fitness_model);
    };

    double compute_phenotype(GenomeStructure const &genome) const {
        double p = genome.background_phenotype;
        for (auto const &locus : loci) { p += genome[locus.first] * (1 + locus.second); }
        return p;
    }

    void parent_loci_sampling(Individual const &parent) {
        for (auto const &locus_pair : parent.loci) {
            if (locus_pair.second or bernouilli_distr(generator)) {
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
        poisson_distribution<u_long> poisson_distr(genome.expected_mutations());
        u_long poisson_draw = poisson_distr(generator);

        if (poisson_draw > 0) {
            uniform_int_distribution<u_long> u_distr(0, genome.number_loci - 1);

            for (u_long mutation{0}; mutation < poisson_draw; mutation++) {
                auto locus = u_distr(generator);
                auto it = loci.find(locus);
                if (it != loci.end()) {
                    if (it->second) {
                        it->second = false;
                    } else {
                        loci.erase(it);
                    }
                    phenotype -= genome.at(locus);
                } else {
                    loci[locus] = false;
                    phenotype += genome.at(locus);
                }
            }
        }
        update_fitness(fitness_model);
    };
};

class Population {
  public:
    // Parameters
    u_long population_size{};
    u_long elapsed{0};
    u_long id{0};
    GenomeStructure genome;
    FitnessModel &fitness_function;

    // Individuals composing the population
    vector<Individual> parents{};
    ~Population() = default;
    explicit Population(u_long population_size, GenomeStructure &genome, FitnessModel &f_model)
        : population_size{population_size}, genome(genome), fitness_function{f_model} {
        parents.resize(population_size);
    };

    bool check() const {
        return all_of(parents.begin(), parents.end(), [this](auto const &p) {
            return abs(p.phenotype - p.compute_phenotype(genome)) < 1e-3;
        });
    }

    void mutation() {
        for (auto &p : parents) { p.mutation(genome, fitness_function); }
    };

    void selection_and_random_mating() {
        vector<double> fitnesses(population_size, 0);
        transform(parents.begin(), parents.end(), fitnesses.begin(),
            [](Individual const &b) { return b.fitness; });
        discrete_distribution<u_long> parent_distr(fitnesses.begin(), fitnesses.end());

        vector<Individual> children;
        children.reserve(parents.size());
        for (u_long i{0}; i < parents.size(); i++) {
            u_long mum = parent_distr(generator);
            u_long dad = parent_distr(generator);
            children.emplace_back(parents[mum], parents[dad], genome, fitness_function);
        }
        parents = move(children);
    };

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
        cout << loci_tmp.size() << " fixation events." << endl;
        for (auto const &fixation : loci_tmp) {
            genome.background_phenotype += genome[fixation] * 2;
            genome[fixation] = -genome[fixation];
            for (auto &p : parents) { p.loci.erase(fixation); }
        }
        assert(check());
    }

    void run(u_long nbr_generations, Trace &trace) {
        cout << "Population " << id << " for " << nbr_generations << " generations." << endl;
        // Run under selection
        u_long last_pct = 0;
        for (u_long sample{1}; sample <= nbr_generations; sample++) {
            elapsed++;
            trace.add("Population", id);
            trace.add("Generation", elapsed);

            TimeVar t_start = timeNow();
            selection_and_random_mating();
            timer.selection += duration(timeNow() - t_start);

            trace.add("Mutational expected var", genome.expected_variance());
            vector<double> phenotypes(population_size, 0);
            transform(parents.begin(), parents.end(), phenotypes.begin(),
                [](Individual const &b) { return b.phenotype; });

            double var_mutational = -variance(phenotypes);
            t_start = timeNow();
            mutation();
            timer.mutation += duration(timeNow() - t_start);

            transform(parents.begin(), parents.end(), phenotypes.begin(),
                [](Individual const &b) { return b.phenotype; });
            var_mutational += variance(phenotypes);
            trace.add("Mutational var", var_mutational);
            trace.add("Phenotype mean", mean(phenotypes));
            trace.add("Phenotype var", variance(phenotypes));

            vector<double> fitnesses(population_size, 0);
            transform(parents.begin(), parents.end(), fitnesses.begin(),
                [](Individual const &b) { return b.fitness; });
            trace.add("Fitness mean", mean(fitnesses));
            trace.add("Fitness var", variance(fitnesses));

            u_long pct = 10 * sample / nbr_generations;
            fitness_function.update();
            if (pct != last_pct) {
                cout << pct * 10 << "% (" << sample << " generations)" << endl;
                last_pct = pct;
                t_start = timeNow();
                fixation();
                timer.fixation += duration(timeNow() - t_start);
            }
        }
    };
};

class Process {
  private:
    vector<Population> populations;

  public:
    explicit Process(
        Population const &pop, u_long nbr_generations, u_long nbr_lineages, Trace &trace) {
        double optimum = pop.fitness_function.get_optimum();
        for (u_long i = 1; i < nbr_lineages + 1; ++i) {
            populations.emplace_back(pop);
            populations.back().id = i;
            populations.back().fitness_function.set_optimum(optimum);
            populations.back().run(nbr_generations, trace);
        }
        timer_cout();
    }

    static void timer_cout() {
        double total_time = timer.mutation + timer.selection + timer.fixation;
        std::cout << std::setprecision(3) << total_time / 1e9 << "s total time" << std::endl;
        std::cout << 100 * timer.mutation / total_time << "% of time spent in mutation ("
                  << timer.mutation / 1e9 << "s)" << std::endl;
        std::cout << 100 * timer.selection / total_time << "% of time spent in selection ("
                  << timer.selection / 1e9 << "s)" << std::endl;
        std::cout << 100 * timer.fixation / total_time << "% of time spent in fixation ("
                  << timer.fixation / 1e9 << "s)" << std::endl;
    }
};