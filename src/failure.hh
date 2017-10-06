#pragma once

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "reliability.hh"
#include "trace.hh"

class FailureMechanism
{
  public:
    const std::string name;

    FailureMechanism(const std::string& _n) : name(_n) {}
    virtual double timeToFailure(const DataPoint& data, double fail=std::numeric_limits<double>::signaling_NaN()) const = 0;
    virtual std::shared_ptr<ReliabilityDistribution> distribution(const std::vector<MTTFSegment>&) const = 0;
};

class NBTI : public FailureMechanism
{
  private:
    const double q = 1.60217662e-19;    // C
    const double k_B = 8.6173303e-5;    // eV/K

    // Device parameters
    const double Vt0 = 0.49158;         // V, taken from PTM [1] at 45nm
    const double Cox = 5.934e-6;        // F
    const double alpha = 1.3;           // alpha power law [2]

    // Low-level parameters (chosen from [3])
    const double A = 9e11;
    const double B = 8.5e11;
    const double C = 3e16;
    const double Gamma_IT = 2.2;
    const double Gamma_HT = 2.2;
    const double Gamma_OT = 9;
    const double E_Akf = 0.175;         // eV
    const double E_Akr = 0.2;           // eV
    const double E_ADH2 = 0.58;         // eV
    const double E_AHT = 0.03;          // eV
    const double E_AOT = 0.15;          // eV
    const double eta = 5e12;
    const double beta_OT = 0.36;

    // High-level parameters
    const double beta = 2;              // Weibull shape parameter [4]
    const double fail_default = 0.03;   // Relative delay change [5]

    const double dt = 3600*24;          // days

  public:
    NBTI() : FailureMechanism("NBTI") {}
    static const std::shared_ptr<FailureMechanism> model() { static NBTI nbti; return std::make_shared<NBTI>(nbti); }
    double degradation(double t, double vdd, double dVth, double temperature, double duty_cycle) const;
    double timeToFailure(const DataPoint& data, double fail=std::numeric_limits<double>::signaling_NaN()) const override;
    std::shared_ptr<ReliabilityDistribution> distribution(const std::vector<MTTFSegment>& mttfs) const override;
};

/*
 * [1] PTM
 * 
 * [2] K. Joshi, S. Mukhopadhyay, N. Goel, and S. Mahapatra, “A consistent
 *     physical framework for N and P BTI in HKMG MOSFETs,” in Reliability
 *     Physics Symposium (IRPS), 2012 IEEE International, 2012, p. 5A.3.1-5A.3.10.
 * 
 * [3]
 * 
 * [4] "Failure Mechanisms and Models for Semiconductor Devices,"" JEDEC Solid
 *     State Technology Institution, JEP122H, Oct. 2011.
 * 
 * [5] F. Oboril and M. B. Tahoori, “ExtraTime: Modeling and analysis of
 *     wearout due to transistor aging at microarchitecture-level,” in
 *     IEEE/IFIP International Conference on Dependable Systems and Networks
 *     (DSN 2012), 2012, pp. 1–12.
 */