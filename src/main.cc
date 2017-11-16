#include <algorithm>
#include <cmath>
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

    int n;
    bool print_rates, v;
    string time_units, dist_file;
    vector<shared_ptr<FailureMechanism>> mechanisms;
    xml_document doc;

    try
    {
        vector<string> units{"seconds", "minutes", "hours", "days", "weeks", "months", "years"};
        ValuesConstraint<string> unit_values(units);

        CmdLine cmd("Compute the reliability distribution of a chip", ' ', "0.1");
        SwitchArg verbose("v", "verbose", "Display progress output", cmd);
        SwitchArg rates("", "print-aging-rates", "Print aging rate of each unit for each trace", cmd);
        ValueArg<string> phenomena("", "aging-mechanisms", "Comma-separated list of aging mechanisms to include or \"all\" for all of them", false, "all", "mechanisms", cmd);
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
        v = verbose.getValue();

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
    }
    catch (ArgException& e)
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }

    if (v)
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
    if (v)
        cout << "Creating failure dependency graph..." << endl;
    shared_ptr<Component> root = make_shared<Group>(doc.child("group"), units);

    if (v)
        cout << "Computing aging rates..." << endl;
    for (const shared_ptr<Unit>& unit: units)
        unit->computeReliability(mechanisms);

    // Monte Carlo sim to get overall failure distribution
    for (int i = 0; i < n; i++)
    {
        if (v)
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