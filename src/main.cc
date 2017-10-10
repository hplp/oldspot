#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <pugixml.hpp>
#include <random>
#include <string>
#include <tclap/CmdLine.h>
#include <utility>
#include <vector>

#include "failure.hh"
#include "trace.hh"
#include "unit.hh"

using namespace std;

int main(int argc, char* argv[])
{
    using namespace pugi;
    using namespace TCLAP;

    int n;
    xml_document doc;

    try
    {
        CmdLine cmd("Compute the reliability distribution of a chip", ' ', "0.1");
        ValueArg<int> iterations("n", "iterations", "Number of Monte-Carlo iterations to perform (default: 1)", false, 1, "iterations", cmd);
        ValueArg<char> delimiter("", "trace-delimiter", "One-character delimiter for data in input trace files (default: ,)", false, ',', "delim", cmd);
        UnlabeledValueArg<string> config("chip-config", "File containing chip configuration", true, "", "filename", cmd);
        cmd.parse(argc, argv);

        xml_parse_result result = doc.load_file(config.getValue().c_str());
        if (!result)
        {
            cerr << config.getValue() << ": " << result.description()
                 << " at offset " << result.offset << endl;
            return 1;
        }

        Unit::delim = delimiter.getValue();
        n = iterations.getValue();
    }
    catch (ArgException& e)
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    vector<shared_ptr<Unit>> units;
    Group root(doc, units);

    vector<shared_ptr<FailureMechanism>> mechanisms = {NBTI::model()};
    for (const shared_ptr<Unit>& unit: units)
        unit->computeReliability(mechanisms);

    // Monte Carlo sim to get overall failure distribution
    vector<double> ttfs(n);
    random_device dev;
    mt19937 gen(dev());
    for (int i = 0; i < n; i++)
    {
        double total_time = 0;
        double fail_time = 0;
        double t_eq_prev = 0;
        for (const shared_ptr<Unit>& unit: units)
        {
            unit->healthy(true);
            unit->current_reliability = 1;
        }
        while (root.healthy())
        {
            // Determine which failure distributions to use based on which units are still healthy
            shared_ptr<Unit> failed = nullptr;
            for (const shared_ptr<Unit>& unit: units)
            {
                if (unit->healthy())
                {
                    uniform_real_distribution<double> r(0, unit->current_reliability);
                    double t = unit->inverse(r(gen));
                    double t_eq = unit->inverse(unit->current_reliability);
                    t -= t_eq;
                    if (failed == nullptr || fail_time > t)
                    {
                        failed = unit;
                        fail_time = t;
                        t_eq_prev = t_eq;
                    }
                }
            }
            total_time += fail_time;
            failed->healthy(false);
            for (shared_ptr<Unit>& unit: units)
                unit->current_reliability = unit->reliability(fail_time + t_eq_prev);
        }
        ttfs[i] = total_time;
    }

    return 0;
}