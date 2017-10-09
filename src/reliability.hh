#pragma once

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

    double reliability(double t) const;
    double inverse(double r) const;
    double mttf() const;

    WeibullDistribution operator*(const WeibullDistribution& other) const;
    WeibullDistribution& operator=(const WeibullDistribution& other);
    WeibullDistribution& operator*=(const WeibullDistribution& other) { return *this = (*this)*other; }
};