#pragma once

#include <cmath>
#include <vector>

namespace oldspot
{

struct MTTFSegment
{
    double duration;
    double mttf;
};

class WeibullDistribution
{
  private:
    double alpha;
    double beta;

  public:
    static WeibullDistribution estimate(const std::vector<double>& ttfs, double beta=2);

    WeibullDistribution(double a, double b) : alpha(a), beta(b) {}
    WeibullDistribution() : WeibullDistribution(1, 1) {}
    WeibullDistribution(const WeibullDistribution& other) : WeibullDistribution(other.alpha, other.beta) {}
    WeibullDistribution(double b, const std::vector<MTTFSegment>& mttfs);

    double reliability(double t) const { return std::exp(-std::pow(t/alpha, beta)); }
    double inverse(double r) const;
    double mttf() const { return alpha*std::tgamma(1/beta + 1); }
    double rate() const { return alpha; }

    double operator()(double t) const { return reliability(t); }
    WeibullDistribution operator*(const WeibullDistribution& other) const;
    WeibullDistribution& operator=(const WeibullDistribution& other);
    WeibullDistribution& operator*=(const WeibullDistribution& other) { return *this = (*this)*other; }
};

} // namespace oldspot