#pragma once

#include <random>

// Random generator engine with seed 0.
u_long seed{0};
std::default_random_engine generator(seed);

std::normal_distribution<double> normal_distrib(0.0, 1.0);
std::bernoulli_distribution bernouilli_distrib(0.5);