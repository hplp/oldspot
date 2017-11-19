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

const shared_ptr<FailureMechanism> NBTI::model()
{
    static shared_ptr<NBTI> nbti = make_shared<NBTI>();
    return nbti;
}

double NBTI::degradation(double t, double vdd, double dVth, double temperature, double duty_cycle) const
{
    duty_cycle = pow(duty_cycle/(1 + sqrt((1 - duty_cycle)/2)), 1.0/6.0);
    double V = vdd - Vt0 - dVth;
    if (V < 0)
    {
        cerr << "warning: subthreshold VDD " << vdd << " not supported" << endl;
        cerr << "         operating at threshold instead" << endl;
        V = 0;
    }
    double E_AIT = 2.0/3.0*(E_Akf - E_Akr) + E_ADH2/6;
    double dN_IT = A*pow(V, Gamma_IT)*exp(-E_AIT/(k_B*temperature))*pow(t, 1.0/6.0);
    double dN_HT = B*pow(V, Gamma_HT)*exp(-E_AHT/(k_B*temperature));

    return duty_cycle*0.027e-12*(dN_IT + dN_HT);
}

double NBTI::timeToFailure(const DataPoint& data, double duty_cycle, double fail) const
{
    if (isnan(fail))
        fail = fail_default;
    if (duty_cycle == 0)
        return numeric_limits<double>::infinity();

    // Create a linear approximation of dVth(t)
    double dVth_fail = (data.data.at("vdd") - Vt0) - (data.data.at("vdd") - Vt0)/pow(1 + fail, 1/alpha); // [ExtraTime]
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

const shared_ptr<FailureMechanism> EM::model()
{
    static shared_ptr<EM> em = make_shared<EM>();
    return em;
}

double EM::timeToFailure(const DataPoint& data, double, double) const
{
    return A*pow(data.data.at("power")/data.data.at("vdd")/(w*h), -n)*exp(Ea/(k_B*data.data.at("temperature")));
}

const shared_ptr<FailureMechanism> HCI::model()
{
    static shared_ptr<HCI> hci = make_shared<HCI>();
    return hci;
}

double HCI::timeToFailure(const DataPoint& data, double duty_cycle, double fail) const
{
    if (isnan(fail))
        fail = fail_default;
    double vdd = data.data.at("vdd");
    double dVth_fail = (vdd - Vt0_n) - (vdd - Vt0_n)/pow(1 + fail, 1/alpha); // [ExtraTime]

    double Vt = k_B/eV_J*data.data.at("temperature")/q;
    double vdsat = ((vdd - Vt0_n + 2*Vt)*L*Esat)/(vdd - Vt0_n + 2*Vt + A_bulk*L*Esat);
    double Em = (vdd - vdsat)/l;
    double Eox = (vdd - Vt0_n)/tox;
    double A_HCI = q/Cox*K*sqrt(Cox*(vdd - Vt0_n));
    double t = pow(dVth_fail/(A_HCI*exp(Eox/E0)*exp(-phi_it/eV_J/(q*lambda*Em))), 1/n)/(duty_cycle*data.data.at("frequency"));

    return t;
}

const shared_ptr<FailureMechanism> TDDB::model()
{
    static shared_ptr<TDDB> tddb = make_shared<TDDB>();
    return tddb;
}

double TDDB::timeToFailure(const DataPoint& data, double, double) const
{
    double T = data.data.at("temperature");
    return pow(data.data.at("vdd"), a - b*T)*exp((X + Y/T + Z*T)/(k_B*T));
}

} // namespace oldspot