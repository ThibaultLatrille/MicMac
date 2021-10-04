#pragma once

#include <numeric>
#include <vector>

// Function for sum
double sum(std::vector<double> const &v) { return std::accumulate(v.begin(), v.end(), 0.0); }

// Function for average
double mean(std::vector<double> const &v) { return sum(v) / static_cast<double>(v.size()); }

// Function for variance
double variance(std::vector<double> const &v) {
    double s2 =
        std::accumulate(v.begin(), v.end(), 0.0, [](double a, double b) { return a + b * b; });
    double m = mean(v);
    return (s2 / static_cast<double>(v.size())) - m * m;
}

// Function for variance
double sum_squarred(std::vector<double> const &v) {
    return std::accumulate(v.begin(), v.end(), 0.0, [](double a, double b) { return a + b * b; });
}