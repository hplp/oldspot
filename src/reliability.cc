#include "reliability.hh"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

WeibullDistribution::WeibullDistribution(double _b, const vector<MTTFSegment>& mttfs)
    : alpha(0), beta(_b)
{
    // Convert MTTFs into rate parameters
    vector<double> alphas(mttfs.size());
    for (size_t i = 0; i < mttfs.size(); i++)
        alphas[i] = mttfs[i].mttf/tgamma(1/beta + 1);
    
    // Accumulate rates into average rate [1]
    for (size_t i = 0; i < alphas.size(); i++)
        alpha += mttfs[i].duration/alphas[i];
    alpha /= accumulate(mttfs.begin(), mttfs.end(), 0,
                        [](double a, MTTFSegment b){ return a + b.duration; });
    
    // Invert to resemble actual Weibull alpha
    alpha = 1/alpha;
}

double WeibullDistribution::reliability(double t) const
{
    return exp(-pow(t/alpha, beta));
}

double WeibullDistribution::inverse(double r) const
{
    return alpha*pow(-log(r), 1/beta);
}

double WeibullDistribution::mttf() const
{
    return alpha*tgamma(1/beta + 1);
}

/*
 * [1] Y. Xiang, T. Chantem, R. P. Dick, X. S. Hu and L. Shang, "System-
 * level reliability modeling for MPSoCs," 2010 IEEE/ACM/IFIP International
 * Conference on Hardware/Software Codesign and System Synthesis (CODES+ISSS),
 * Scottsdale, AZ, 2010, pp. 297-306.
 */