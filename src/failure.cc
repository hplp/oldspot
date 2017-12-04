#include "failure.hh"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "reliability.hh"
#include "trace.hh"
#include "util.hh"

namespace oldspot
{

using namespace std;

FailureMechanism::FailureMechanism(const string& _n) : name(_n)
{
    p = {{"L", 65},         // nm
         {"Vt0_p", 0.5},    // PMOS threshold voltage, V
         {"Vt0_n", 0.5},    // NMOS threshold voltage, V
         {"tox", 1.8},      // nm
         {"Cox", 1.92e-20}, // F/nm^2
         {"alpha", 1.3}};   // alpha power law [2]
}

NBTI::NBTI() : FailureMechanism("NBTI")
{
    p["A"] = 5.5e12;
    p["B"] = 8e11;
    p["Gamma_IT"] = 4.5;
    p["Gamma_HT"] = 4.5;
    p["E_Akf"] = 0.175; // eV
    p["E_Akr"] = 0.2;   // eV
    p["E_ADH2"] = 0.58; // eV
    p["E_AHT"] = 0.03;  // eV
}

double NBTI::degradation(double t, double vdd, double dVth, double temperature, double duty_cycle) const
{
    duty_cycle = pow(duty_cycle/(1 + sqrt((1 - duty_cycle)/2)), 1.0/6.0);
    double V = vdd - p.at("Vt0_p") - dVth;
    if (V < 0)
    {
        cerr << "warning: subthreshold VDD " << vdd << " not supported" << endl;
        cerr << "         operating at threshold instead" << endl;
        V = 0;
    }
    double E_AIT = 2.0/3.0*(p.at("E_Akf") - p.at("E_Akr")) + p.at("E_ADH2")/6;
    double dN_IT = p.at("A")*pow(V, p.at("Gamma_IT"))*exp(-E_AIT/(k_B*temperature))*pow(t, 1.0/6.0);
    double dN_HT = p.at("B")*pow(V, p.at("Gamma_HT"))*exp(-p.at("E_AHT")/(k_B*temperature));

    return duty_cycle*0.027e-12*(dN_IT + dN_HT);
}

double NBTI::timeToFailure(const DataPoint& data, double duty_cycle, double fail) const
{
    if (isnan(fail))
        fail = fail_default;
    if (duty_cycle == 0)
        return numeric_limits<double>::infinity();

    // Create a linear approximation of dVth(t)
    double dVth_fail = (data.data.at("vdd") - p.at("Vt0_p")) - (data.data.at("vdd") - p.at("Vt0_p"))/pow(1 + fail, 1/p.at("alpha")); // [ExtraTime]
    double dVth = 0, dVth_prev = 0;
    double t = 0;
    for (; dVth < dVth_fail; t += dt)
    {
        dVth_prev = dVth;
        dVth = degradation(t, data.data.at("vdd"), dVth, data.data.at("temperature"), duty_cycle);
    }
    t -= dt;

    if (dVth == 0)
        return 0;
    else // Linearly interpolate to find time at which dVth == dVth_fail (MTTF)
        return linterp(dVth_fail, {dVth_prev, t - dt}, {dVth, t});
}

EM::EM() : FailureMechanism("EM")
{
    p["n"] = 2;
    p["Ea"] = 0.8;      // eV
    p["w"] = 4.5e-7;    // m
    p["h"] = 1.2e-6;    // m
    p["A"] = 3.22e21;   // extracted from [?]
}

double EM::timeToFailure(const DataPoint& data, double, double) const
{
    return p.at("A")*pow(data.data.at("power")/data.data.at("vdd")/(p.at("w")*p.at("h")), -p.at("n"))*exp(p.at("Ea")/(k_B*data.data.at("temperature")));
}

HCI::HCI() : FailureMechanism("HCI")
{
    p["E0"] = 0.8;      // V/nm
    p["K"] = 1.7e8;     // nm/C^0.5
    p["A_bulk"] = 0.005;
    p["phi_it"] = 3.7;  // eV
    p["lambda"] = 7.8;  // nm
    p["l"] = 17;        // nm
    p["Esat"] = 0.011;  // V/nm
    p["n"] = 0.45;
}

double HCI::timeToFailure(const DataPoint& data, double duty_cycle, double fail) const
{
    if (isnan(fail))
        fail = fail_default;
    double vdd = data.data.at("vdd");
    double dVth_fail = (vdd - p.at("Vt0_n")) - (vdd - p.at("Vt0_n"))/pow(1 + fail, 1/p.at("alpha")); // [ExtraTime]

    double Vt = k_B/eV_J*data.data.at("temperature")/q;
    double vdsat = ((vdd - p.at("Vt0_n") + 2*Vt)*p.at("L")*p.at("Esat"))/(vdd - p.at("Vt0_n") + 2*Vt + p.at("A_bulk")*p.at("L")*p.at("Esat"));
    double Em = (vdd - vdsat)/l;
    double Eox = (vdd - p.at("Vt0_n"))/p.at("tox");
    double A_HCI = q/p.at("Cox")*p.at("K")*sqrt(p.at("Cox")*(vdd - p.at("Vt0_n")));
    double t = pow(dVth_fail/(A_HCI*exp(Eox/p.at("E0"))*exp(-p.at("phi_it")/eV_J/(q*p.at("lambda")*Em))), 1/p.at("n"))/(duty_cycle*data.data.at("frequency"));

    return t;
}

double TDDB::timeToFailure(const DataPoint& data, double, double) const
{
    double T = data.data.at("temperature");
    return pow(data.data.at("vdd"), a - b*T)*exp((X + Y/T + Z*T)/(k_B*T));
}

} // namespace oldspot

/*
 * [2] K. Joshi, S. Mukhopadhyay, N. Goel, and S. Mahapatra, “A consistent
 *     physical framework for N and P BTI in HKMG MOSFETs,” in Reliability
 *     Physics Symposium (IRPS), 2012 IEEE International, 2012, p. 5A.3.1-5A.3.10.
 */