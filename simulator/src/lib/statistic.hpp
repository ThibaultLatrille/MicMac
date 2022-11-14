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
// Slope of a linear regression
double slope(const std::vector<double> &x, const std::vector<double> &y) {
    assert(x.size() == y.size());
    const auto n = static_cast<double>(x.size());
    const double s_x = std::accumulate(x.begin(), x.end(), 0.0);
    const double s_y = std::accumulate(y.begin(), y.end(), 0.0);
    const double s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const double s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const double a = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    return a;
}

class Stat {
  public:
    double n{0};
    double total{0};

    Stat() = default;
    ~Stat() = default;

    void add(double const &x, double const &w = 1.0) {
        total += x * w;
        n += w;
    }

    double mean() const { return total / n; }
};
