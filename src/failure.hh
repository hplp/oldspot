#pragma once

#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "reliability.hh"
#include "trace.hh"

namespace oldspot
{

class FailureMechanism
{
  protected:
    const double beta = 2; // Weibull shape parameter [4]

  public:
    // Universal constants
    static constexpr double q = 1.60217662e-19;     // C
    static constexpr double k_B = 8.6173303e-5;     // eV/K
    static constexpr double eV_J = 6.242e18;        // eV -> J
    static constexpr double nm2_cm2 = 1e14;         // nm^2 -> cm^2
    // Device parameters
    static constexpr double L = 65;                 // nm
    static constexpr double Vt0 = 0.5;            // V
    static constexpr double Vt0_n = 0.5;          // V
    static constexpr double tox = 1.8;              // nm
    static constexpr double Cox = 1.92e-20;         // F/nm^2
    static constexpr double alpha = 1.3;            // alpha power law [2]

    const double fail_default = 0.05;   // Relative delay change [5]

    const std::string name;

    FailureMechanism(const std::string& _n) : name(_n) {}
    virtual double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const = 0;
    virtual WeibullDistribution distribution(const std::vector<MTTFSegment>& mttfs) const
    {
        return WeibullDistribution(beta, mttfs);
    }
};

class NBTI : public FailureMechanism
{
  private:
    // Low-level parameters (chosen from [3])
    const double A = 5.5e12;
    const double B = 8e11;
    const double Gamma_IT = 4.5;
    const double Gamma_HT = 4.5;
    const double E_Akf = 0.175;         // eV
    const double E_Akr = 0.2;           // eV
    const double E_ADH2 = 0.58;         // eV
    const double E_AHT = 0.03;          // eV

    // High-level parameters 
    const double dt = 3600*24;          // days

  public:
    NBTI() : FailureMechanism("NBTI") {}
    double degradation(double t, double vdd, double dVth, double temperature, double duty_cycle) const;
    double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

class EM : public FailureMechanism
{
  private:
    // Low-level parameters (chosen from [7])
    const double n = 2;
    const double Ea = 0.8;      // eV
    const double w = 4.5e-7;    // m
    const double h = 1.2e-6;    // m
    const double A = 3.22e21;   // Extracted from [6]

  public:
    EM() : FailureMechanism("EM") {}
    double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

class HCI : public FailureMechanism
{
  private:
    // Low-level parameters (chosen from [8])
    const double E0 = 0.8;              // V/nm
    const double K = 1.7e8;             // nm/C^0.5
    const double A_bulk = 0.005;
    const double phi_it = 3.7;          // eV
    const double lambda = 7.8;          // nm
    const double l = 17;                // nm
    const double Esat = 0.011;          // V/nm
    const double n = 0.45;

  public:
    HCI() : FailureMechanism("HCI") {}
    double degradation(double t, double vdd, double temperature, double frequency, double duty_cycle) const;
    double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

class TDDB : public FailureMechanism
{
  private:
    // Low-level parameters (chosen from [9])
    const double a = 78;
    const double b = -0.081;    // 1/K
    const double X = 0.759;     // eV
    const double Y = -66.8;     // eV*K
    const double Z = -8.37e-4;  // eV/K

  public:
    TDDB() : FailureMechanism("TDDB") {}
    double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

} // namespace oldspot

namespace std
{
    template<>
    struct less<shared_ptr<oldspot::FailureMechanism>>
    {
        bool operator()(const shared_ptr<oldspot::FailureMechanism>& lhs, const shared_ptr<oldspot::FailureMechanism>& rhs) const
        {
            return lhs->name < rhs->name;
        }
    };
}

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
 * 
 * [6]
 * 
 * [7]
 * 
 * [8]
 * 
 * [9]
 */