#include "failure.hh"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "interp.hh"
#include "reliability.hh"
#include "trace.hh"

namespace oldspot
{

using namespace std;

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
    double n = eta*pow(V, -Gamma_OT/beta_OT)*exp(E_AOT/(k_B*temperature*beta_OT));
    double dN_OT = C*(1 - exp(-pow(t/n, beta_OT)));
    return duty_cycle*q/Cox*(dN_IT + dN_HT + dN_OT);
}

double NBTI::timeToFailure(const DataPoint& data, double fail) const
{
    if (isnan(fail))
        fail = fail_default;

    // Create a linear approximation of dVth(t)
    double dVth_fail = (data.data.at("vdd") - Vt0) - (data.data.at("vdd") - Vt0)/pow(1 + fail, 1/alpha); // [ExtraTime]
    double dVth = 0, dVth_prev = 0;
    double t = 0;
    for (; dVth < dVth_fail; t += dt)
    {
        dVth_prev = dVth;
        dVth = degradation(t, data.data.at("vdd"), dVth, data.data.at("temperature"), data.data.at("activity"));
    }
    t -= dt;

    if (dVth == 0)
        return 0;
    else // Linearly interpolate to find time at which dVth == dVth_fail (MTTF)
        return linterp(dVth_fail, {dVth_prev, t - dt}, {dVth, t});
}

double EM::timeToFailure(const DataPoint& data, double fail) const
{
    return A*pow(data.data.at("power")/data.data.at("vdd")/(w*h), -n)*exp(Ea/(k_B*data.data.at("temperature")));
}

} // namespace oldspot