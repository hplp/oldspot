#include "failure.hh"

#include <cmath>
#include <fstream>
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

/**
 * Constructor for all FailureMechanisms.  Initializes process-dependent parameters
 * that are used across all failure mechanisms.  Default parameter values come from:
 * [1] R. Vattikonda, W. Wang, Y. Cao, "Modeling and minimization of PMOS NBTI effect
 *     for robust nanometer design," DAC, pp. 1047-1052, 2006.
 */
FailureMechanism::FailureMechanism(const string& _n, const string& tech_file) : name(_n)
{
    // Default values
    p = {{"L", 65},         // nm
         {"Vt0_p", 0.5},    // PMOS threshold voltage, V
         {"Vt0_n", 0.5},    // NMOS threshold voltage, V
         {"tox", 1.8},      // nm
         {"Cox", 1.92e-20}, // F/nm^2
         {"alpha", 1.3}};   // alpha power law [2]
    if (!tech_file.empty())
    {
        Parameters params = read_params(tech_file);
        p.insert(params.begin(), params.end());
    }
}

/**
 * Read parameters from a file.  The file should consist of name-value pairs
 * separated by tabs, with one pair on each line.  Lines beginning with '#'
 * are comments and ignored by the parser.
 */
FailureMechanism::Parameters
FailureMechanism::read_params(const string& file)
{
    Parameters params;
    fstream f(file);
    if (f)
    {
        string line;
        while (getline(f, line))
        {
            if (line[0] != '#')
            {
                vector<string> tokens = split(line, '\t');
                if (tokens.size() != 2)
                    warn("%s: %d: unable to parse line\n", file.c_str(), line);
                else
                    params[tokens[0]] = stod(tokens[1]);
            }
        }
    }
    else
        warn("%s: file not found\n", file.c_str());
    return params;
}

/**
 * Constructor for the NBTI aging mechanism.  Adds parameters to the parameter list
 * that are specific to NBTI, which are read from a parameter file (see read_params).
 * Default parameters come from:
 * [2] Joshi, K., Mukhopadhyay, S., Goel, N., and Mahapatra, S. "A consistent physical
 *     framework for N and P BTI in HKMG MOSFETs." In Reliability Physics Symposium (IRPS),
 *     2012 IEEE International, pages 5A{3. IEEE, 2012.
 */
NBTI::NBTI(const string& tech_file, const string& nbti_file) : FailureMechanism("NBTI", tech_file)
{
    // Default values
    p["A"] = 5.5e12;
    p["B"] = 8e11;
    p["Gamma_IT"] = 4.5;
    p["Gamma_HT"] = 4.5;
    p["E_Akf"] = 0.175; // eV
    p["E_Akr"] = 0.2;   // eV
    p["E_ADH2"] = 0.58; // eV
    p["E_AHT"] = 0.03;  // eV
    if (!nbti_file.empty())
    {
        Parameters params = read_params(nbti_file);
        p.insert(params.begin(), params.end());
    }
}

/**
 * Compute the degradation in threshold voltage due to NBTI over a given period
 * of time at the given temperature and duty cycle.  The model for this degradation
 * comes from [2] (see NBTI::NBTI()).
 */
double
NBTI::degradation(double t, double vdd, double dVth, double temperature, double duty_cycle) const
{
    duty_cycle = pow(duty_cycle/(1 + sqrt((1 - duty_cycle)/2)), 1.0/6.0);
    double V = vdd - p.at("Vt0_p") - dVth;
    if (V < 0)
    {
        warn("subthreshold VDD %f not supported; operating at threshold instead\n", vdd);
        V = 0;
    }
    double E_AIT = 2.0/3.0*(p.at("E_Akf") - p.at("E_Akr")) + p.at("E_ADH2")/6;
    double dN_IT = p.at("A")*pow(V, p.at("Gamma_IT"))*exp(-E_AIT/(k_B*temperature))*pow(t, 1.0/6.0);
    double dN_HT = p.at("B")*pow(V, p.at("Gamma_HT"))*exp(-p.at("E_AHT")/(k_B*temperature));

    return duty_cycle*0.027e-12*(dN_IT + dN_HT);
}

/**
 * Estimate the time to failure for NBTI.  Since the model from [2] is not invertible,
 * this is done by finding a start and end time that contains the point of failure and
 * then linearly interpolating between them.  The amount of time each segment takes is
 * defined by NBTI::dt.  The effective duty cycle for NBTI is computed from:
 * [3] F. Oboril and M. B. Tahoori, "ExtraTime: Modeling and analysis of
 *     wearout due to transistor aging at microarchitecture-level," in
 *     IEEE/IFIP International Conference on Dependable Systems and Networks
 *     (DSN 2012), 2012, pp. 1–12.
 */
double
NBTI::timeToFailure(const DataPoint& data, double duty_cycle, double fail) const
{
    if (isnan(fail))
        fail = fail_default;
    if (duty_cycle == 0)
        return numeric_limits<double>::infinity();

    // Create a linear approximation of dVth(t)
    double dVth_fail = (data.data.at("vdd") - p.at("Vt0_p"))
                       - (data.data.at("vdd") - p.at("Vt0_p"))/pow(1 + fail, 1/p.at("alpha")); // [3]
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

/**
 * Constructor for the EM aging mechanism.  Adds parameters to the parameter list
 * that are specific to EM, which are read from a parameter file (see read_params).
 * Default parameters come from:
 * [4] 
 */
EM::EM(const string& tech_file, const string& em_file) : FailureMechanism("EM", tech_file)
{
    p["n"] = 2;
    p["Ea"] = 0.8;          // eV
    p["w"] = 4.5e-7;        // m
    p["h"] = 1.2e-6;        // m
    p["A"] = 3.22e21;       // extracted from [?]
    p["wire_density"] = 1;  // wires/m^2
    if (!em_file.empty())
    {
        Parameters params = read_params(em_file);
        p.insert(params.begin(), params.end());
    }
}

/**
 * Compute mean-time-to-failure of EM using Black's Equation:
 * [5] Black, J. R.  Electromigration--a brief survey and some recent results.
 *     IEEE Transactions on Electron Devices, 16(4):338–347, 1969.
 */
double
EM::timeToFailure(const DataPoint& data, double, double) const
{
    double j = 0;
    if (data.data.count("current_density") != 0)
        j = data.data.at("current_density");
    else if (data.data.count("current") != 0)
        j = data.data.at("current")/(p.at("w")*p.at("h"));
    else
    {
        warn("current density or current not found in trace data; approximating as P/V\n");
        j = data.data.at("power")/data.data.at("vdd")/(p.at("w")*p.at("h"));
    }
    return p.at("A")*pow(j, -p.at("n"))*exp(p.at("Ea")/(k_B*data.data.at("temperature")));
}

/**
 * Constructor for the HCI aging mechanism.  Adds parameters to the parameter list
 * that are specific to HCI, which are read from a parameter file (see read_params).
 * Default parameters come from [1].
 */
HCI::HCI(const string& tech_file, const std::string& hci_file) : FailureMechanism("HCI", tech_file)
{
    p["E0"] = 0.8;      // V/nm
    p["K"] = 1.7e8;     // nm/C^0.5
    p["A_bulk"] = 0.005;
    p["phi_it"] = 3.7;  // eV
    p["lambda"] = 7.8;  // nm
    p["l"] = 17;        // nm
    p["Esat"] = 0.011;  // V/nm
    p["n"] = 0.45;
    if (!hci_file.empty())
    {
        Parameters params = read_params(hci_file);
        p.insert(params.begin(), params.end());
    }
}

/**
 * Compute the time to failure due to HCI.  The model in [1] is invertible, so
 * just compute the time it takes to reach a failure state.
 */
double
HCI::timeToFailure(const DataPoint& data, double duty_cycle, double fail) const
{
    if (isnan(fail))
        fail = fail_default;
    double vdd = data.data.at("vdd");
    double dVth_fail = (vdd - p.at("Vt0_n")) - (vdd - p.at("Vt0_n"))/pow(1 + fail, 1/p.at("alpha")); // [3]

    double Vt = k_B/eV_J*data.data.at("temperature")/q;
    double vdsat = ((vdd - p.at("Vt0_n") + 2*Vt)*p.at("L")*p.at("Esat"))
                   /(vdd - p.at("Vt0_n") + 2*Vt + p.at("A_bulk")*p.at("L")*p.at("Esat"));
    double Em = (vdd - vdsat)/p.at("l");
    double Eox = (vdd - p.at("Vt0_n"))/p.at("tox");
    double A_HCI = q/p.at("Cox")*p.at("K")*sqrt(p.at("Cox")*(vdd - p.at("Vt0_n")));
    double t = pow(dVth_fail/(A_HCI*exp(Eox/p.at("E0"))*exp(-p.at("phi_it")/eV_J/(q*p.at("lambda")*Em))), 1/p.at("n"))
               /(duty_cycle*data.data.at("frequency"));

    return t;
}

/**
 * Constructor for the TDDB aging mechanism.  Adds parameters to the parameter list
 * that are specific to TDDB, which are read from a parameter file (see read_params).
 * Default parameters come from:
 * [6] Srinivasan, J., Adve, S. V., Bose, P., and Rivers, J. A. "The case for lifetime
 *     reliability-aware microprocessors." In ACM SIGARCH Computer Architecture News,
 *     volume 32, page 276. IEEE Computer Society, 2004.
 */
TDDB::TDDB(const string& tech_file, const string& tddb_file) : FailureMechanism("TDDB", tech_file)
{
    p["a"] = 78;
    p["b"] = -0.081;    // 1/K
    p["X"] = 0.759;     // eV
    p["Y"] = -66.8;     // eV*K
    p["Z"] = -8.37e-4;  // eV/K
    if (!tddb_file.empty())
    {
        Parameters params = read_params(tddb_file);
        p.insert(params.begin(), params.end());
    }
}

/**
 * Compute mean time to failure for TDDB using the equation from [6].
 */
double
TDDB::timeToFailure(const DataPoint& data, double, double) const
{
    double T = data.data.at("temperature");
    return pow(data.data.at("vdd"), p.at("b")*T - p.at("a"))*exp((p.at("X") + p.at("Y")/T + p.at("Z")*T)/(k_B*T));
}

} // namespace oldspot
