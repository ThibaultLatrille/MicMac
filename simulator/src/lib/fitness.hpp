#pragma once

#include "random.hpp"
#include "tclap/CmdLine.h"

class FitnessModel {
  protected:
    double optimum{0};

  public:
    explicit FitnessModel() = default;
    virtual ~FitnessModel() = default;

    virtual double fitness(double v) const = 0;

    void set_optimum(double v) { optimum = v; }
    double get_optimum() const { return optimum; }
    virtual void update() = 0;
};

class NeutralModel : public FitnessModel {
  public:
    explicit NeutralModel() = default;
    ~NeutralModel() override = default;

    double fitness(double v) const override { return 1.0; }
    void update() override{};
};

class GaussianModel : public FitnessModel {
  protected:
    double peakness{};
    double epistasis{};

  public:
    explicit GaussianModel() = default;
    explicit GaussianModel(double const &peakness, double const &epistasis)
        : peakness{peakness}, epistasis{epistasis} {}

    double fitness(double v) const override { return exp(-peakness * pow(v, epistasis)); }
    void update() override{};
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

    GaussianModel get_model(u_long nb_loci) { return GaussianModel(peakness.getValue() / static_cast<double>(nb_loci), epistasis.getValue()); }
};

class MovingOptimumModel : public GaussianModel {
  protected:
    double std_optimum{0};

  public:
    explicit MovingOptimumModel() = default;
    explicit MovingOptimumModel(
        double const &peakness, double const &epistasis, double const &std_optimum)
        : GaussianModel(peakness, epistasis), std_optimum{std_optimum} {}

    double fitness(double v) const override {
        return exp(-peakness * pow((v - optimum), epistasis));
    }
    void update() override { optimum += std_optimum * normal_distrib(generator); };

    ~MovingOptimumModel() override = default;
};

class MovingOptimumArgParse : GaussianArgParse {
  public:
    explicit MovingOptimumArgParse(TCLAP::CmdLine &cmd) : GaussianArgParse(cmd) {}

    TCLAP::ValueArg<double> std_optimum{"", "std_optimum",
        "Standard deviation of the change in optimum.", false, 1.0, "double", cmd};

    MovingOptimumModel get_moving_optimum_model(u_long nb_loci) {
        return MovingOptimumModel(
            peakness.getValue() / static_cast<double>(nb_loci), epistasis.getValue(), std_optimum.getValue());
    }
};

class DirectionalModel final : public FitnessModel {
  public:
    explicit DirectionalModel() = default;
    ~DirectionalModel() override = default;

    double fitness(double v) const override {
        if (v < 0) {
            return 0.0;
        } else {
            return v;
        }
    }
    void update() override{};
};
