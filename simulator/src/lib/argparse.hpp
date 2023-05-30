#pragma once

#include <cassert>
#include <cmath>
#include <fstream>
#include "io.hpp"
#include "tclap/CmdLine.h"

class OutputArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit OutputArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}
    TCLAP::ValueArg<double> number_of_generations{"", "number_of_generations",
        "Number of generations from root to leaf", false, 1000, "u_long", cmd};
    TCLAP::ValueArg<u_long> number_of_lineages{
        "", "number_of_lineages", "Number of lineages", false, 2, "u_long", cmd};
    TCLAP::ValueArg<std::string> tree_path{
        "", "tree", "input tree path (newick or nhx format)", false, "", "string", cmd};
    TCLAP::ValueArg<u_long> burn_in{"", "burn_in",
        "Burn-in period in number of generations (before speciation).", false, 1000, "u_long", cmd};
    TCLAP::ValueArg<std::string> output_path{"o", "output", "output path", true, "", "string", cmd};
    TCLAP::SwitchArg ouput_tsv{
        "", "ouput_tsv", "Output also as a tsv file (per generation)", cmd, false};
    TCLAP::ValueArg<u_long> seed{
        "", "seed", "Random number generator seed", false, 0, "u_long", cmd};
    TCLAP::ValueArg<u_long> seed_pop_size{"", "seed_pop_size",
        "Random number generator seed (specific to the population size)", false, 0, "u_long", cmd};
    TCLAP::ValueArg<u_long> seed_mut_rate{"", "seed_mut_rate",
        "Random number generator seed (specific to the mutation rate)", false, 0, "u_long", cmd};
};
