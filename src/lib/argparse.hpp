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

    TCLAP::ValueArg<u_long> population_size{
        "", "population_size", "Number of individuals", false, 100, "u_long", cmd};
    TCLAP::ValueArg<u_long> number_of_generations{"", "number_of_generations",
        "Number of year between generations", false, 1000, "u_long", cmd};
    TCLAP::ValueArg<std::string> output_path{"o", "output", "output path", true, "", "string", cmd};
    TCLAP::ValueArg<u_long> seed{
        "", "seed", "Random number generation seed", false, 0, "u_long", cmd};
};
