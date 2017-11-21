#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <pugixml.hpp>
#include <random>
#include <string>
#include <tclap/CmdLine.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "failure.hh"
#include "trace.hh"
#include "unit.hh"
#include "util.hh"

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

    vector<string> time_units{"seconds", "minutes", "hours", "days", "weeks", "months", "years"};
    ValuesConstraint<string> time_constraint(time_units);

    vector<shared_ptr<FailureMechanism>> mechanisms;
    xml_document doc;

    CmdLine cmd("Compute the reliability distribution of a chip", ' ', "0.1");
    SwitchArg verbose("v", "verbose", "Display progress output", cmd);
    SwitchArg separate("", "separate-aging-rates", "Display a second table with aging rates separated by mechanism per unit (only works for fresh configuration)", cmd);
    ValueArg<string> phenomena("", "aging-mechanisms", "Comma-separated list of aging mechanisms to include or \"all\" for all of them", false, "all", "mechanisms", cmd);
    ValueArg<char> delimiter("", "trace-delimiter", "One-character delimiter for data in input trace files (default: ,)", false, ',', "delim", cmd);
    ValueArg<string> time("", "time-units", "Units for displaying time to failure (default: hours)", false, "hours", &time_constraint, cmd);
    ValueArg<int> iterations("n", "iterations", "Number of Monte-Carlo iterations to perform (default: 1000)", false, 1000, "iterations", cmd);
    ValueArg<string> dist_dump("", "dump-ttfs", "Dump time-to-failure distribution to file", false, "", "filename", cmd);
    UnlabeledValueArg<string> config("chip-config", "File containing chip configuration", true, "", "filename", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (ArgException& e)
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    xml_parse_result result = doc.load_file(config.getValue().c_str());
    if (!result)
    {
        cerr << config.getValue() << ": " << result.description()
                << " at offset " << result.offset << endl;
        return 1;
    }

    Unit::delim = delimiter.getValue();

    string p = phenomena.getValue();
    transform(p.begin(), p.end(), p.begin(), ::tolower);
    if (p == "all")
        mechanisms = {NBTI::model(), EM::model(), HCI::model(), TDDB::model()};
    else
    {
        for (const string& token: split(p, ','))
        {
            if (token == "nbti")
                mechanisms.push_back(NBTI::model());
            else if (token == "em")
                mechanisms.push_back(EM::model());
            else if (token == "hci")
                mechanisms.push_back(HCI::model());
            else if (token == "tddb")
                mechanisms.push_back(TDDB::model());
            else
                cerr << "warning: ignoring unknown aging mechanism \"" << token << '"' << endl;
        }
    }
    if (mechanisms.empty())
    {
        cerr << "error: no aging mechanisms selected" << endl;
        return 1;
    }

    if (verbose.getValue())
        cout << "Creating units..." << endl;
    vector<shared_ptr<Unit>> units;
    for (const xml_node& child: doc.children("unit"))
    {
        if (node_is(child, "unit"))
            units.push_back(make_shared<Unit>(child, units.size()));
        else if (node_is(child, "core"))
            units.push_back(make_shared<Core>(child, units.size()));
        else if (node_is(child, "logic"))
            units.push_back(make_shared<Logic>(child, units.size()));
        else if (node_is(child, "memory"))
            units.push_back(make_shared<Memory>(child, units.size()));
        else
        {
            cerr << "unknown unit type \"" << child.name() << '"' << endl;
            exit(1);
        }
    }
    if (verbose.getValue())
        cout << "Creating failure dependency graph..." << endl;
    shared_ptr<Component> root = make_shared<Group>(doc.child("group"), units);

    if (verbose.getValue())
        cout << "Computing aging rates..." << endl;
    for (const shared_ptr<Unit>& unit: units)
        unit->compute_reliability(mechanisms);

    // Monte Carlo sim to get overall failure distribution
    for (int i = 0; i < iterations.getValue(); i++)
    {
        if (verbose.getValue())
            cout << "Beginning Monte Carlo iteration " << i << endl;

        unordered_set<shared_ptr<Component>> failed_components;
        unordered_set<shared_ptr<Unit>> healthy(units.begin(), units.end());
        double t = 0;
        for (shared_ptr<Unit>& unit: units)
            unit->reset();
        while (!root->failed())
        {
            for (const shared_ptr<Unit>& unit: units)
                if (!unit->failed())
                    unit->set_configuration(root);

            double dt_event = numeric_limits<double>::infinity();
            shared_ptr<Unit> failed;
            for (const shared_ptr<Unit>& unit: healthy)
            {
                double dt = unit->get_next_event();
                if (dt_event > dt)
                {
                    failed = unit;
                    dt_event = dt;
                }
            }
            
            if (isinf(dt_event))
                break;

            for (const shared_ptr<Unit>& unit: healthy)
                unit->update_reliability(dt_event);
            failed->failure();
            if (failed->failed())
                healthy.erase(healthy.find(failed));
            t += dt_event;

            Component::walk(root, [&](const shared_ptr<Component>& c) {
                if (c->failed() && failed_components.count(c) == 0)
                {
                    c->ttfs.push_back(t);
                    failed_components.insert(c);
                }
            });
            for (const shared_ptr<Unit>& unit: Unit::parents_failed(root, units))
                failed_components.insert(unit);
        }
    }

    sort(units.begin(), units.end(), [](const shared_ptr<Unit>& a, const shared_ptr<Unit>& b) {
        return b->ttfs.size() < a->ttfs.size();
    });
    vector<string> rows{root->name};
    vector<string> cols{"MTTF", "Failures", "Aging Rate"};
    transform(units.begin(), units.end(), back_inserter(rows), [](const shared_ptr<Unit>& u){ return u->name; });
    unordered_map<string, unordered_map<string, double>> data;
    data[root->name]["MTTF"] = convert_time(root->mttf(), time.getValue());
    data[root->name]["Failures"] = root->ttfs.size();
    data[root->name]["Aging Rate"] = numeric_limits<double>::quiet_NaN();
    for (const shared_ptr<Unit>& u: units)
    {
        data[u->name]["MTTF"] = convert_time(u->mttf(), time.getValue());
        data[u->name]["Failures"] = u->ttfs.size();
        data[u->name]["Aging Rate"] = convert_time(u->aging_rate(Unit::fresh), time.getValue());
    }
    print_table(rows, cols, data);

    if (separate.getValue())
    {
        sort(units.begin(), units.end(), [&](const shared_ptr<Unit>& a, const shared_ptr<Unit>& b) {
            return b->aging_rate(Unit::fresh) < a->aging_rate(Unit::fresh);
        });
        vector<string> rows;
        vector<string> cols;
        transform(units.begin(), units.end(), back_inserter(rows), [](const shared_ptr<Unit>& u){ return u->name; });
        transform(mechanisms.begin(), mechanisms.end(), back_inserter(cols), [](const shared_ptr<FailureMechanism>& m){ return m->name; });
        data.clear();
        for (const shared_ptr<Unit>& u: units)
            for (const shared_ptr<FailureMechanism>& m: mechanisms)
                data[u->name][m->name] = convert_time(u->aging_rate(m), time.getValue());
        cout << endl;
        print_table(rows, cols, data);
    }

    if (!dist_dump.getValue().empty())
    {
        ofstream dist(dist_dump.getValue());
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