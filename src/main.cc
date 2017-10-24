#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
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
#include <unordered_set>
#include <utility>
#include <vector>

#include "failure.hh"
#include "trace.hh"
#include "unit.hh"

using namespace oldspot;
using namespace pugi;
using namespace std;

inline bool node_is(const xml_node& node, const string& type)
{
    return strcmp(node.attribute("type").value(), type.c_str()) == 0;
}

// Assumes time is already in seconds
inline double convert_time(double time, const string& units)
{
    if (units == "seconds")
        return time;
    time /= 60;
    if (units == "minutes")
        return time;
    time /= 60;
    if (units == "hours")
        return time;
    time /= 24;
    if (units == "days")
        return time;
    time /= 7;
    if (units == "weeks")
        return time;
    time /= 4;
    if (units == "months")
        return time;
    time /= 12;
    if (units == "years")
        return time;
    throw invalid_argument("unknown time unit \"" + units + '"');
}

int main(int argc, char* argv[])
{
    using namespace TCLAP;

    int n;
    bool print_rates;
    string time_units, dist_file;
    xml_document doc;

    try
    {
        vector<string> units = {"seconds", "minutes", "hours", "days", "weeks", "months", "years"};
        ValuesConstraint<string> unit_values(units);

        CmdLine cmd("Compute the reliability distribution of a chip", ' ', "0.1");
        SwitchArg rates("", "print-aging-rates", "Print aging rate of each unit for each trace", cmd);
        ValueArg<char> delimiter("", "trace-delimiter", "One-character delimiter for data in input trace files (default: ,)", false, ',', "delim", cmd);
        ValueArg<string> time("", "time-units", "Units for displaying time to failure (default: hours)", false, "hours", &unit_values, cmd);
        ValueArg<int> iterations("n", "iterations", "Number of Monte-Carlo iterations to perform (default: 1000)", false, 1000, "iterations", cmd);
        ValueArg<string> dist_dump("", "dump-ttfs", "Dump time-to-failure distribution to file", false, "", "filename", cmd);
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
        time_units = time.getValue();
        print_rates = rates.getValue();
        dist_file = dist_dump.getValue();
    }
    catch (ArgException& e)
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    vector<shared_ptr<Unit>> units;
    for (const xml_node& child: doc.children("unit"))
    {
        if (node_is(child, "unit"))
            units.push_back(make_shared<Unit>(child, units.size()));
        else if (node_is(child, "core"))
            units.push_back(make_shared<Core>(child, units.size()));
        else if (node_is(child, "logic"))
            units.push_back(make_shared<Logic>(child, units.size()));
        else
        {
            cerr << "unknown unit type \"" << child.name() << '"' << endl;
            exit(1);
        }
    }
    shared_ptr<Component> root = make_shared<Group>(doc.child("group"), units);

    vector<shared_ptr<FailureMechanism>> mechanisms = {NBTI::model()};
    for (const shared_ptr<Unit>& unit: units)
        unit->computeReliability(mechanisms);

    Unit::init_configurations(root, units);
    if (print_rates)
    {
        const string f = "(failed)";
        size_t rate_width = f.length();
        size_t name_width = 0;
        for (const shared_ptr<Unit>& unit: units)
        {
            name_width = max(name_width, unit->name.length());
            for (size_t i = 0; i < Unit::configurations(); i++)
                rate_width = max(rate_width, to_string(convert_time(unit->aging_rate(i), time_units)).length());
        }
        cout << "Aging Rates:" << endl;
        for (const shared_ptr<Unit>& unit: units)
        {
            cout << left << setw(name_width) << unit->name << " | ";
            for (size_t i = 0; i < Unit::configurations(); i++)
            {
                cout << right << setw(rate_width) << (unit->aging_rate(i) == 0 ? f : to_string(convert_time(unit->aging_rate(i), time_units)));
                if (i != Unit::configurations() - 1)
                    cout << " | ";
            }
            cout << endl;
        }
        cout << endl;
    }

    // Monte Carlo sim to get overall failure distribution
    random_device dev;
    mt19937 gen(dev());
    for (int i = 0; i < n; i++)
    {
        unordered_set<shared_ptr<Component>> failed_components;
        double total_time = 0;
        double fail_time = 0;
        double t_eq_prev = 0;
        for (shared_ptr<Unit>& unit: units)
            unit->reset();
        while (!root->failed())
        {
            shared_ptr<Unit> failed = nullptr;
            Unit::set_configuration(units);
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
            failed->add_failure();
            for (const shared_ptr<Unit>& unit: units)
                unit->current_reliability = unit->reliability(fail_time + t_eq_prev);
            if (failed->serial())
                failed->current_reliability = 1;

            Component::walk(root, [&](shared_ptr<Component>& c) {
                if (c->failed() && failed_components.count(c) == 0)
                {
                    c->ttfs.push_back(total_time);
                    failed_components.insert(c);
                }
            });
        }
    }

    cout << "MTTFs:" << endl;
    cout << root->name << ": " << convert_time(root->mttf(), time_units)
                       << " (" << convert_time(root->mttf_interval().first, time_units)
                       << ',' << convert_time(root->mttf_interval().second, time_units) << ')' << endl;
    for (const shared_ptr<Unit>& unit: units)
        cout << unit->name << ": " << convert_time(unit->mttf(), time_units)
                           << " (" << convert_time(unit->mttf_interval().first, time_units)
                           << ',' << convert_time(unit->mttf_interval().second, time_units) << ')' << endl;

    if (!dist_file.empty())
    {
        ofstream dist(dist_file);
        if (dist)
        {
            dist << root->name;
            for (double ttf: root->ttfs)
                dist << ',' << ttf;
            dist << endl;
            for (const shared_ptr<Unit>& unit: units)
            {
                dist << unit->name;
                for (double ttf: unit->ttfs)
                    dist << ',' << ttf;
                dist << endl;
            }
        }
    }

    return 0;
}