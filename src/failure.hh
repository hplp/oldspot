#pragma once

#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "reliability.hh"
#include "trace.hh"

namespace oldspot
{

class FailureMechanism
{
  protected:
    const double beta = 2; // Weibull shape parameter [4]
    std::unordered_map<std::string, double> p; // Device parameters

    std::unordered_map<std::string, double> read_params(const std::string& file);

  public:
    // Universal constants
    static constexpr double q = 1.60217662e-19;     // C
    static constexpr double k_B = 8.6173303e-5;     // eV/K
    static constexpr double eV_J = 6.242e18;        // eV -> J
    static constexpr double nm2_cm2 = 1e14;         // nm^2 -> cm^2

    static constexpr double fail_default = 0.05;    // Relative delay change [5]

    const std::string name;

    FailureMechanism(const std::string& _n, const std::string& tech_file);
    virtual double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const = 0;
    virtual WeibullDistribution distribution(const std::vector<MTTFSegment>& mttfs) const
    {
        return WeibullDistribution(beta, mttfs);
    }
};

class NBTI : public FailureMechanism
{
  private:
    static constexpr double dt = 3600*24;          // days

  public:
    NBTI(const std::string& tech_file, const std::string& nbti_file);
    double degradation(double t, double vdd, double dVth, double temperature, double duty_cycle) const;
    double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

class EM : public FailureMechanism
{
  public:
    EM(const std::string& tech_file, const std::string& em_file);
    double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

class HCI : public FailureMechanism
{
  public:
    HCI(const std::string& tech_file, const std::string& hci_file);
    double degradation(double t, double vdd, double temperature, double frequency, double duty_cycle) const;
    double timeToFailure(const DataPoint& data, double duty_cycle, double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

class TDDB : public FailureMechanism
{
  public:
    TDDB(const std::string& tech_file, const std::string& tddb_file);
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