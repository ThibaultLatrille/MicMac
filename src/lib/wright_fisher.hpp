#include <map>
#include <set>
#include "argparse.hpp"
#include "fitness.hpp"
#include "random.hpp"
#include "statistic.hpp"
using namespace TCLAP;
using namespace std;

class GenomeStructure : public vector<double> {
  public:
    u_long number_loci{};
    double mutation_mean_effect_size{};
    double mutation_rate_per_loci_per_generation{};

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
        phenotype = 0.0;
        for (auto const &locus : loci) { phenotype += genome[locus.first] * (1 + locus.second); }
        update_fitness(fitness_model);
    };

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
    GenomeStructure const &genome{};
    FitnessModel const &fitness_function;

    // Individuals composing the population
    vector<Individual> parents{};
    ~Population() = default;
    explicit Population(
        u_long population_size, GenomeStructure const &genome, FitnessModel const &f_model)
        : population_size{population_size}, genome(genome), fitness_function{f_model} {
        parents.resize(population_size);
    };

    void mutation() {
        for (auto &Individual : parents) { Individual.mutation(genome, fitness_function); }
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
            u_long dad = mum;
            while (dad == mum) { dad = parent_distr(generator); }
            children.emplace_back(parents[mum], parents[dad], genome, fitness_function);
        }
        parents = move(children);
    };

    void run(u_long nbr_generations, Trace &trace) {
        cout << "Population " << id << " for " << nbr_generations << " generations." << endl;
        // Run under selection
        u_long last_pct = 0;
        for (u_long sample{1}; sample <= nbr_generations; sample++) {
            mutation();
            selection_and_random_mating();
            elapsed++;
            add_to_trace(trace);
            u_long pct = 10 * sample / nbr_generations;
            if (pct != last_pct) {
                cout << pct * 10 << "% (" << sample << " generations)" << endl;
                last_pct = pct;
            }
        }
    };

    void add_to_trace(Trace &trace) {
        vector<double> fitnesses(population_size, 0);
        transform(parents.begin(), parents.end(), fitnesses.begin(),
            [](Individual const &b) { return b.fitness; });

        vector<double> phenotypes(population_size, 0);
        transform(parents.begin(), parents.end(), phenotypes.begin(),
            [](Individual const &b) { return b.phenotype; });

        trace.add("Population", id);
        trace.add("Generation", elapsed);
        trace.add("Fitness mean", mean(fitnesses));
        trace.add("Fitness var", variance(fitnesses));
        trace.add("Phenotype mean", mean(phenotypes));
        trace.add("Phenotype var", variance(phenotypes));
    };
};

class Process {
  private:
    Population pop1;
    Population pop2;

  public:
    explicit Process(Population &pop, u_long nbr_generations, Trace &trace) : pop1{pop}, pop2{pop} {
        pop1.id = 1;
        pop1.run(nbr_generations, trace);
        pop2.id = 2;
        pop2.run(nbr_generations, trace);
    }
};