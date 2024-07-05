### **`simulator`**
This folder contains the simulation program to generate the data for the simulated analysis, required to run the experiments in the `data_simulated` folder.

Compiling the simulation program requires the buildsystem `make` and `cmake` and a c++ compiler (`g++` or `clang`).\
Compiling is then done by running the following commands in the `simulator` folder:
```bash
make release
```

This will create the executable `neutral` and `selection` in the `build` folder.\
They can be run with the following command to see the usage:
```bash
build/neutral --help
build/selection --help
```

Usage of the `neutral` program:
```
build/neutral  [--mutation_rate_brownian_sigma <double>]
                  [--mutation_rate_max <double>] [--mutation_rate_min
                  <double>] [--mutation_rate_per_loci_per_generation
                  <double>] [--population_size_noise_theta <double>]
                  [--population_size_noise_sigma <double>]
                  [--population_size_brownian_sigma <double>]
                  [--population_size_max <u_long>] [--population_size_min
                  <u_long>] [--population_size <u_long>]
                  [--variance_environment <double>]
                  [--mutation_mean_effect_size <double>] [--number_loci
                  <u_long>] [--seed_mut_rate <u_long>] [--seed_pop_size
                  <u_long>] [--seed <u_long>] [--ouput_tsv] -o <string>
                  [--burn_in <u_long>] [--tree <string>]
                  [--number_of_lineages <u_long>] [--number_of_generations
                  <u_long>] [--] [--version] [-h]

Where: 

   --mutation_rate_brownian_sigma <double>
     The Brownian sigma (0<sigma) applied to μ at each generation

   --mutation_rate_max <double>
     Maximum mutation_rate

   --mutation_rate_min <double>
     Minimum mutation_rate

   --mutation_rate_per_loci_per_generation <double>
     Mutation rate per locus per generation

   --population_size_noise_theta <double>
     The Ornstein–Uhlenbeck theta (0<=theta<1) applied to Ne at each
     generation

   --population_size_noise_sigma <double>
     The Ornstein–Uhlenbeck sigma (0<sigma) applied to Ne at each
     generation

   --population_size_brownian_sigma <double>
     The Brownian sigma (0<sigma) applied to Ne at each generation

   --population_size_max <u_long>
     Maximum population size

   --population_size_min <u_long>
     Minimum population size

   --population_size <u_long>
     Number of individuals

   --variance_environment <double>
     Environmental variance (Ve).

   --mutation_mean_effect_size <double>
     Mutation mean effect size

   --number_loci <u_long>
     Number of loci

   --seed_mut_rate <u_long>
     Random number generator seed (specific to the mutation rate)

   --seed_pop_size <u_long>
     Random number generator seed (specific to the population size)

   --seed <u_long>
     Random number generator seed

   --ouput_tsv
     Output also as a tsv file (per generation)

   -o <string>,  --output <string>
     (required)  output path

   --burn_in <u_long>
     Burn-in period in number of generations (before speciation).

   --tree <string>
     input tree path (newick or nhx format)

   --number_of_lineages <u_long>
     Number of lineages

   --number_of_generations <u_long>
     Number of generations from root to leaf

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```


Usage of the `selection` program:
```
build/selection  [--mutation_rate_brownian_sigma <double>]
                    [--mutation_rate_max <double>] [--mutation_rate_min
                    <double>] [--mutation_rate_per_loci_per_generation
                    <double>] [--epistasis <double>] [--peakness <double>]
                    [--theta_optimum <double>] [--sigma_optimum <double>]
                    [--bias_optimum <double>]
                    [--population_size_noise_theta <double>]
                    [--population_size_noise_sigma <double>]
                    [--population_size_brownian_sigma <double>]
                    [--population_size_max <u_long>] [--population_size_min
                    <u_long>] [--population_size <u_long>]
                    [--variance_environment <double>]
                    [--mutation_mean_effect_size <double>] [--number_loci
                    <u_long>] [--seed_mut_rate <u_long>] [--seed_pop_size
                    <u_long>] [--seed <u_long>] [--ouput_tsv] -o <string>
                    [--burn_in <u_long>] [--tree <string>]
                    [--number_of_lineages <u_long>]
                    [--number_of_generations <u_long>] [--] [--version]
                    [-h]

Where: 

   --mutation_rate_brownian_sigma <double>
     The Brownian sigma (0<sigma) applied to μ at each generation

   --mutation_rate_max <double>
     Maximum mutation_rate

   --mutation_rate_min <double>
     Minimum mutation_rate

   --mutation_rate_per_loci_per_generation <double>
     Mutation rate per locus per generation

   --epistasis <double>
     'beta' parameter (epistasis) of fitness function
     (exp(-alpha*(phenotype^beta))

   --peakness <double>
     'alpha' parameter (peakness) of the fitness function
     (exp(-alpha*(phenotype^beta))

   --theta_optimum <double>
     The Ornstein–Uhlenbeck theta (0<=theta<1) applied to the optimum at
     each generation

   --sigma_optimum <double>
     The Ornstein–Uhlenbeck sigma (>0) applied to the fitness optimum at
     each generation

   --bias_optimum <double>
     The Ornstein–Uhlenbeck bias (>0) applied to the fitness optimum at
     each generation

   --population_size_noise_theta <double>
     The Ornstein–Uhlenbeck theta (0<=theta<1) applied to Ne at each
     generation

   --population_size_noise_sigma <double>
     The Ornstein–Uhlenbeck sigma (0<sigma) applied to Ne at each
     generation

   --population_size_brownian_sigma <double>
     The Brownian sigma (0<sigma) applied to Ne at each generation

   --population_size_max <u_long>
     Maximum population size

   --population_size_min <u_long>
     Minimum population size

   --population_size <u_long>
     Number of individuals

   --variance_environment <double>
     Environmental variance (Ve).

   --mutation_mean_effect_size <double>
     Mutation mean effect size

   --number_loci <u_long>
     Number of loci

   --seed_mut_rate <u_long>
     Random number generator seed (specific to the mutation rate)

   --seed_pop_size <u_long>
     Random number generator seed (specific to the population size)

   --seed <u_long>
     Random number generator seed

   --ouput_tsv
     Output also as a tsv file (per generation)

   -o <string>,  --output <string>
     (required)  output path

   --burn_in <u_long>
     Burn-in period in number of generations (before speciation).

   --tree <string>
     input tree path (newick or nhx format)

   --number_of_lineages <u_long>
     Number of lineages

   --number_of_generations <u_long>
     Number of generations from root to leaf

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```