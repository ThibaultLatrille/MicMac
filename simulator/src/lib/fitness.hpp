#pragma once

#include "random.hpp"
#include "tclap/CmdLine.h"

class FitnessState {
  public:
    double drift{0};
    double fitness_optimum{0};
    explicit FitnessState() = default;
    ~FitnessState() = default;
};

class FitnessModel {
  protected:
    FitnessState state{};

  public:
    explicit FitnessModel() = default;
    virtual ~FitnessModel() = default;

    virtual double fitness(double v) const = 0;

    void set_state(FitnessState const &s) { state = s; }
    FitnessState get_state() const { return state; }
    virtual void update() = 0;
};

class NeutralModel : public FitnessModel {
  public:
    explicit NeutralModel() = default;
    ~NeutralModel() override = default;

    double fitness(double v) const override { return 1.0; }
    void update() override{};
};


class OrnsteinUhlenbeckBiasModel : public FitnessModel {
  protected:
    double peakness{};
    double epistasis{};
    double bias{};
    double sigma{};
    double theta{};

  public:
    explicit OrnsteinUhlenbeckBiasModel() = default;
    explicit OrnsteinUhlenbeckBiasModel(double const &peakness, double const &epistasis,
        double const &bias_optimum, double const &sigma_optimum, double const &theta_optimum)
        : peakness{peakness},
          epistasis{epistasis},
          bias{bias_optimum},
          sigma{sigma_optimum},
          theta{theta_optimum} {}

    void update() override {
        state.drift += bias;
        state.fitness_optimum +=
            sigma * normal_distrib(generator) + theta * (state.drift - state.fitness_optimum);
    }

    double fitness(double v) const override {
        return exp(-peakness * pow((v - state.fitness_optimum), epistasis));
    }
    ~OrnsteinUhlenbeckBiasModel() override = default;
};

class OrnsteinUhlenbeckBiasArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit OrnsteinUhlenbeckBiasArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}

    TCLAP::ValueArg<double> bias_optimum{"", "bias_optimum",
        "The Ornstein–Uhlenbeck bias (>0) applied to the fitness optimum at each generation", false,
        0.1, "double", cmd};
    TCLAP::ValueArg<double> sigma_optimum{"", "sigma_optimum",
        "The Ornstein–Uhlenbeck sigma (>0) applied to the fitness optimum at each generation",
        false, 0.1, "double", cmd};
    TCLAP::ValueArg<double> theta_optimum{"", "theta_optimum",
        "The Ornstein–Uhlenbeck theta (0<=theta<1) applied to the optimum at each generation",
        false, 0.9, "double", cmd};
    TCLAP::ValueArg<double> peakness{"", "peakness",
        "'alpha' parameter (peakness) of the fitness function "
        "(exp(-alpha*(phenotype^beta))",
        false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> epistasis{"", "epistasis",
        "'beta' parameter (epistasis) of fitness "
        "function (exp(-alpha*(phenotype^beta))",
        false, 2.0, "double", cmd};


    OrnsteinUhlenbeckBiasModel get_model() {
        return OrnsteinUhlenbeckBiasModel(peakness.getValue(), epistasis.getValue(),
            bias_optimum.getValue(), sigma_optimum.getValue(), theta_optimum.getValue());
    }
};
