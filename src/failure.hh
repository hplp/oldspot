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

/**
 * This is the base class for all failure mechanisms.  It defines universal constants
 * and process parameters that apply to all aging mechanisms and ensures that each
 * aging mechanism has the same interface for computing its time to failure.
 * Currently all aging mechanisms are assumed to follow a Weibull distribution with
 * shape parameter 2.
 * 
 * References:
 * [1] "Failure Mechanisms and Models for Semiconductor Devices,"" JEDEC Solid
 *     State Technology Institution, JEP122H, Oct. 2011.
 * [2] F. Oboril and M. B. Tahoori, "ExtraTime: Modeling and analysis of
 *     wearout due to transistor aging at microarchitecture-level," in
 *     IEEE/IFIP International Conference on Dependable Systems and Networks
 *     (DSN 2012), 2012, pp. 1–12.
 */
class FailureMechanism
{
  protected:
    typedef std::unordered_map<std::string, double> Parameters;

    const double beta = 2; // Weibull shape parameter [1]
    Parameters p; // Device parameters

    Parameters read_params(const std::string& file);

  public:
    // Universal constants
    static constexpr double q = 1.60217662e-19;     // C
    static constexpr double k_B = 8.6173303e-5;     // eV/K
    static constexpr double eV_J = 6.242e18;        // Convert eV -> J

    static constexpr double fail_default = 0.05;    // Relative delay change [2]

    const std::string name;

    FailureMechanism(const std::string& _n, const std::string& tech_file);
    virtual double timeToFailure(const DataPoint& data, double duty_cycle,
                                 double fail=std::numeric_limits<double>::signaling_NaN()) const = 0;
    virtual WeibullDistribution distribution(const std::vector<MTTFSegment>& mttfs) const
    {
        return WeibullDistribution(beta, mttfs);
    }
};

/**
 * Definition for failures due to negative bias temperature instability (NBTI).
 */
class NBTI : public FailureMechanism
{
  private:
    static constexpr double dt = 3600*24; // days

  public:
    NBTI(const std::string& tech_file, const std::string& nbti_file);
    double degradation(double t, double vdd, double dVth, double temperature, double duty_cycle) const;
    double timeToFailure(const DataPoint& data, double duty_cycle,
                         double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

/**
 * Definition for failures due to electromigration (EM).
 */
class EM : public FailureMechanism
{
  public:
    EM(const std::string& tech_file, const std::string& em_file);
    double timeToFailure(const DataPoint& data, double duty_cycle,
                         double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

/**
 * Definition for failures due to hot-carrier injection (HCI).
 */
class HCI : public FailureMechanism
{
  public:
    HCI(const std::string& tech_file, const std::string& hci_file);
    double degradation(double t, double vdd, double temperature, double frequency, double duty_cycle) const;
    double timeToFailure(const DataPoint& data, double duty_cycle,
                         double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

/**
 * Definition for failures due to time-dependent dielectric breakdown (TDDB).
 */
class TDDB : public FailureMechanism
{
  public:
    TDDB(const std::string& tech_file, const std::string& tddb_file);
    double timeToFailure(const DataPoint& data, double duty_cycle,
                         double fail=std::numeric_limits<double>::signaling_NaN()) const override;
};

} // namespace oldspot

namespace std
{

/**
 * Comparison functor to allow FailureMechanisms to be stored in ordered
 * std::maps and std::sets and be sorted by name.
 */
template<>
struct less<shared_ptr<oldspot::FailureMechanism>>
{
    bool operator()(const shared_ptr<oldspot::FailureMechanism>& lhs,
                    const shared_ptr<oldspot::FailureMechanism>& rhs) const
    {
        return lhs->name < rhs->name;
    }
};

}