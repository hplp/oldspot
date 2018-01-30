#pragma once

#include <cmath>
#include <vector>

namespace oldspot
{

/**
 * Pair of values representing a time period over which a device experiences a failure
 * rate with a particular mean time to failure.
 */
struct MTTFSegment
{
    double duration;
    double mttf;
};

/**
 * The Weibull distribution is a method for representing the failure probability of a
 * device over time (or, equivalently, the the fraction of surviving devices within a
 * population) where the failure rate increases with time, as it does with most aging
 * mechanisms.  A Weibull distribution takes the form of:
 * 
 *    R(t) = exp(-(t/a)^b)
 * 
 * where a is the rate parameter and b is the shape parameter (most often represented
 * as alpha and beta).  Currently, aging mechanisms are all assumed to follow a
 * Weibull distribution with beta = 2 and alpha dependent on mechanism and operating
 * conditions (voltage, temperature, duty cycle, etc.) [1]. This class also has
 * helper methods for computing the Weibull distribution of a sytem with components
 * that have their own Weibull distributions or when the rate parameter changes
 * over time.
 * 
 * References:
 * [1] "Failure Mechanisms and Models for Semiconductor Devices,"" JEDEC Solid
 *     State Technology Institution, JEP122H, Oct. 2011.
 */
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