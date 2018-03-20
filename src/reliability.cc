#include "reliability.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace oldspot
{

using namespace std;

/**
 * Estimate the rate parameter of the Weibull-distributed times-to-failure with
 * the given shape parameter.
 */
WeibullDistribution
WeibullDistribution::estimate(const std::vector<double>& ttfs, double beta)
{
    WeibullDistribution dist;
    dist.beta = beta;
    dist.alpha = pow(accumulate(ttfs.begin(), ttfs.end(), 0.0, [&](double a, double b){
        return a + b*b;
    })/ttfs.size(), 1/beta);
    return dist;
}

/**
 * Create a Weibull distribution using the given set of time-varying mean-times-to-failure,
 * computed using:
 * [1] Y. Xiang, T. Chantem, R. P. Dick, X. S. Hu and L. Shang, "System-
 *     level reliability modeling for MPSoCs," 2010 IEEE/ACM/IFIP International
 *     Conference on Hardware/Software Codesign and System Synthesis (CODES+ISSS),
 *     Scottsdale, AZ, 2010, pp. 297-306.
 */
WeibullDistribution::WeibullDistribution(double b, const vector<MTTFSegment>& mttfs)
    : WeibullDistribution(1, b)
{
    // Convert MTTFs into rate parameters
    vector<double> alphas;
    alphas.reserve(mttfs.size());
    transform(mttfs.begin(), mttfs.end(), back_inserter(alphas),
              [this](const MTTFSegment& a){ return a.mttf/tgamma(1/beta + 1); });
    
    // Accumulate rates into average rate [1]
    alpha = 0.0;
    double total_time = 0.0;
    for (const MTTFSegment& mttf: mttfs)
    {
        alpha += mttf.duration/mttf.mttf;
        total_time += mttf.duration;
    }
    alpha /= total_time;

    // Invert to resemble actual Weibull alpha
    alpha = 1/alpha;
}

/**
 * Compute the time it takes to get to a particular reliabilty value with this
 * Weibull distribution's parameters.
 */
double
WeibullDistribution::inverse(double r) const
{
    if (isinf(alpha))
        return numeric_limits<double>::infinity();
    return alpha*pow(-log(r), 1/beta);
}

/**
 * Compute the resulting Weibull distribution when multiplying two Weibull
 * distributions with the same shape parameter (if they have different shape
 * parameters, the result doesn't follow the Weibull distribution)
 */
WeibullDistribution
WeibullDistribution::operator*(const WeibullDistribution& other) const
{
    if (beta != other.beta)
        throw invalid_argument("the product of two Weibull distributions with different shapes does not follow a Weibull distribution");
    double a = pow(pow(1/alpha, beta) + pow(1/other.alpha, beta), -1/beta);
    return WeibullDistribution(a, beta);
}

/**
 * Copy the parameters of the given Weibull distribution.
 */
WeibullDistribution&
WeibullDistribution::operator=(const WeibullDistribution& other)
{
    alpha = other.alpha;
    beta = other.beta;
    return *this;
}

} // namespace oldspot