#include "reliability.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace oldspot
{

using namespace std;

WeibullDistribution WeibullDistribution::estimate(const std::vector<double>& ttfs, double beta)
{
    WeibullDistribution dist;
    dist.beta = beta;
    dist.alpha = accumulate(ttfs.begin(), ttfs.end(), 0.0, [&](double a, double b){
        return a + b*b;
    })/ttfs.size();
    return dist;
}

WeibullDistribution::WeibullDistribution(double b, const vector<MTTFSegment>& mttfs)
    : WeibullDistribution(1, b)
{
    // Convert MTTFs into rate parameters
    vector<double> alphas;
    alphas.reserve(mttfs.size());
    transform(mttfs.begin(), mttfs.end(), back_inserter(alphas),
              [this](const MTTFSegment& a){ return a.mttf/tgamma(1/beta + 1); });
    
    // Accumulate rates into average rate [1]
    alpha = inner_product(mttfs.begin(), mttfs.end(), alphas.begin(), 0.0,
                          plus<double>(),
                          [](const MTTFSegment& m, double a){ return m.duration/a; })/
            accumulate(mttfs.begin(), mttfs.end(), 0.0,
                       [](double a, MTTFSegment b){ return a + b.duration; });

    // Invert to resemble actual Weibull alpha
    alpha = 1/alpha;
}

double WeibullDistribution::inverse(double r) const
{
    if (isinf(alpha))
        return numeric_limits<double>::infinity();
    return alpha*pow(-log(r), 1/beta);
}

WeibullDistribution WeibullDistribution::operator*(const WeibullDistribution& other) const
{
    if (beta != other.beta)
        throw invalid_argument("the product of two Weibull distributions with different shapes does not follow a Weibull distribution");
    double a = pow(pow(1/alpha, beta) + pow(1/other.alpha, beta), -1/beta);
    return WeibullDistribution(a, beta);
}

WeibullDistribution& WeibullDistribution::operator=(const WeibullDistribution& other)
{
    alpha = other.alpha;
    beta = other.beta;
    return *this;
}

} // namespace oldspot

/*
 * [1] Y. Xiang, T. Chantem, R. P. Dick, X. S. Hu and L. Shang, "System-
 * level reliability modeling for MPSoCs," 2010 IEEE/ACM/IFIP International
 * Conference on Hardware/Software Codesign and System Synthesis (CODES+ISSS),
 * Scottsdale, AZ, 2010, pp. 297-306.
 */