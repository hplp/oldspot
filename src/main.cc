#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
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
    bool print_rates;
    xml_document doc;

    try
    {
        CmdLine cmd("Compute the reliability distribution of a chip", ' ', "0.1");
        SwitchArg rates("", "print-aging-rates", "Print aging rate of each unit for each trace", cmd);
        ValueArg<char> delimiter("", "trace-delimiter", "One-character delimiter for data in input trace files (default: ,)", false, ',', "delim", cmd);
        ValueArg<int> iterations("n", "iterations", "Number of Monte-Carlo iterations to perform (default: 1000)", false, 1000, "iterations", cmd);
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
        print_rates = rates.getValue();
    }
    catch (ArgException& e)
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    vector<shared_ptr<Unit>> units;
    shared_ptr<Component> root = make_shared<Group>(doc, units, n);

    vector<shared_ptr<FailureMechanism>> mechanisms = {NBTI::model()};
    for (const shared_ptr<Unit>& unit: units)
        unit->computeReliability(mechanisms);

    for (uint64_t i = 0; i < (1ULL << units.size()); i++)
    {
        for (const shared_ptr<Unit>& unit: units)
            unit->failed((i&(1 << unit->id)) != 0);
        if (!root->failed())
            Unit::add_configuration(i);
    }
    if (print_rates)
    {
        const string f = "(failed)";
        size_t rate_width = f.length();
        size_t name_width = 0;
        for (const shared_ptr<Unit>& unit: units)
        {
            name_width = max(name_width, unit->name.length());
            for (size_t i = 0; i < Unit::configurations(); i++)
                rate_width = max(rate_width, to_string(unit->aging_rate(i)).length());
        }
        for (const shared_ptr<Unit>& unit: units)
        {
            cout << left << setw(name_width) << unit->name << " | ";
            for (size_t i = 0; i < Unit::configurations(); i++)
            {
                cout << right << setw(rate_width) << (unit->aging_rate(i) == 0 ? f : to_string(unit->aging_rate(i)));
                if (i != Unit::configurations() - 1)
                    cout << " | ";
            }
            cout << endl;
        }
    }

    // Monte Carlo sim to get overall failure distribution
    random_device dev;
    mt19937 gen(dev());
    for (int i = 0; i < n; i++)
    {
        double total_time = 0;
        double fail_time = 0;
        double t_eq_prev = 0;
        for (shared_ptr<Unit>& unit: units)
            unit->reset();
        Component::walk(root, [&](const shared_ptr<Component>& c){
            c->ttfs[i] = numeric_limits<double>::infinity();
        });
        while (!root->failed())
        {
            Unit::set_configuration(units);
            shared_ptr<Unit> failed = nullptr;
            for (const shared_ptr<Unit>& unit: units)
            {
                if (!unit->failed())
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
            failed->failed(true);
            for (const shared_ptr<Unit>& unit: units)
                unit->current_reliability = unit->reliability(fail_time + t_eq_prev);
            
            Component::walk(root, [&](const shared_ptr<Component>& c) {
                if (c->failed())
                    c->ttfs[i] = min(c->ttfs[i], total_time);
            });
        }
        Component::walk(root, [&](const shared_ptr<Component>& c) {
            if (isinf(c->ttfs[i]))
                c->ttfs[i] = root->ttfs[i];
        });
    }

    Component::walk(root, [](const shared_ptr<Component>& c) {
        cout << c->name << ": " << c->mttf()/(60*60*24*365) << endl;
    });

    return 0;
}