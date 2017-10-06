#pragma once

#include <vector>

struct MTTFSegment
{
    double duration;
    double mttf;
};

class ReliabilityDistribution
{
  public:
    virtual double reliability(double t) const = 0; // aka survivor function
    virtual double inverse(double r) const = 0; // time to reach reliability r
    virtual double mttf() const = 0;
    // pdf, CDF, hazard, etc.
};

class WeibullDistribution : public ReliabilityDistribution
{
  private:
    double alpha;
    double beta;

  public:
    WeibullDistribution(double _b, const std::vector<MTTFSegment>& mttfs);

    double reliability(double t) const override;
    double inverse(double r) const override;
    double mttf() const override;
};