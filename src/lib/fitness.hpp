#pragma once

#include "tclap/CmdLine.h"

class FitnessModel {
  public:
    explicit FitnessModel() = default;

    virtual double fitness(double v) const = 0;

    virtual ~FitnessModel() = default;
};

class NeutralModel : public FitnessModel {
  public:
    explicit NeutralModel() = default;

    double fitness(double v) const override { return 1.0; }

    ~NeutralModel() override = default;
};

class GaussianModel final : public FitnessModel {
  private:
    double peakness{};
    double epistasis{};

  public:
    explicit GaussianModel() = default;
    explicit GaussianModel(double const &peakness, double const &epistasis)
        : peakness{peakness}, epistasis{epistasis} {}

    double fitness(double v) const override { return exp(-peakness * pow(v, epistasis)); }

    ~GaussianModel() override = default;
};

class GaussianArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit GaussianArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}

    TCLAP::ValueArg<double> peakness{"", "peakness",
        "'alpha' parameter (peakness) of the fitness function "
        "(exp(-alpha*(phenotype^beta))",
        false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> epistasis{"", "epistasis",
        "'beta' parameter (epistasis) of fitness "
        "function (exp(-alpha*(phenotype^beta))",
        false, 2.0, "double", cmd};

    GaussianModel get_model() { return GaussianModel(peakness.getValue(), epistasis.getValue()); }
};