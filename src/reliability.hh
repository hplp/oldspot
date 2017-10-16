#pragma once

#include <cmath>
#include <vector>

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
    WeibullDistribution(double a, double b) : alpha(a), beta(b) {}
    WeibullDistribution() : WeibullDistribution(1, 1) {}
    WeibullDistribution(const WeibullDistribution& other) : WeibullDistribution(other.alpha, other.beta) {}
    WeibullDistribution(double b, const std::vector<MTTFSegment>& mttfs);

    double reliability(double t) const { return std::exp(-std::pow(t/alpha, beta)); }
    double inverse(double r) const { return alpha*std::pow(-std::log(r), 1/beta); }
    double mttf() const { return alpha*std::tgamma(1/beta + 1); }
    double rate() const { return alpha; }

    WeibullDistribution operator*(const WeibullDistribution& other) const;
    WeibullDistribution& operator=(const WeibullDistribution& other);
    WeibullDistribution& operator*=(const WeibullDistribution& other) { return *this = (*this)*other; }
};