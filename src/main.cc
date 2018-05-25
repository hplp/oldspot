#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <pugixml.hpp>
#include <random>
#include <set>
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

inline bool
node_is(const xml_node& node, const string& type)
{
    return strcmp(node.attribute("type").value(), type.c_str()) == 0;
}

// Assumes time is already in seconds
inline double
convert_time(double time, const string& units)
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

int
main(int argc, char* argv[])
{
    using namespace TCLAP;

    vector<string> time_units{"seconds", "minutes", "hours", "days", "weeks", "months", "years"};
    ValuesConstraint<string> time_constraint(time_units);

    set<shared_ptr<FailureMechanism>> mechanisms;
    xml_document doc;

    CmdLine cmd("Compute the reliability distribution of a chip", ' ', "0.1");
    SwitchArg verbose("v", "verbose", "Display progress output", cmd);
    ValueArg<string> tddb("", "tddb-parameters", "File containing model parameters for TDDB", false, "", "filename", cmd);
    ValueArg<string> hci("", "hci-parameters", "File containing model parameters for HCI", false, "", "filename", cmd);
    ValueArg<string> em("", "em-parameters", "File containing model parameters for electromigration", false, "", "filename", cmd);
    ValueArg<string> nbti("", "nbti-parameters", "File containing model parameters for NBTI", false, "", "filename", cmd);
    ValueArg<string> technology("", "technology-file", "File containing technology constants for aging mechanisms", false, "", "filename", cmd);
    ValueArg<string> phenomena("", "aging-mechanisms", "Comma-separated list of aging mechanisms to include or \"all\" for all of them", false, "all", "mechanisms", cmd);
    ValueArg<char> delimiter("", "trace-delimiter", "One-character delimiter for data in input trace files (default: ,)", false, ',', "delim", cmd);
    ValueArg<string> time("", "time-units", "Units for displaying time to failure (default: hours)", false, "hours", &time_constraint, cmd);
    ValueArg<string> separate("", "mechanism-aging-rates", "Write per-mechanism aging rates for each unit to file (only works for fresh configuration)", false, "", "filename", cmd);
    ValueArg<string> dist_dump("", "dump-ttfs", "Dump time-to-failure distribution to file", false, "", "filename", cmd);
    ValueArg<string> rates("", "unit-aging-rates", "Write per-unit aging rates, MTTFs, and failure counts to file (aging rates only for fresh configuration)", false, "", "filename", cmd);
    ValueArg<int> iterations("n", "iterations", "Number of Monte-Carlo iterations to perform (default: 1000)", false, 1000, "iterations", cmd);
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

    transform(phenomena.getValue().begin(), phenomena.getValue().end(), phenomena.getValue().begin(), ::tolower);
    for (const string& token: split(phenomena.getValue(), ','))
    {
        if (token == "nbti" || token == "all")
            mechanisms.insert(make_shared<NBTI>(technology.getValue(), nbti.getValue()));
        if (token == "em" || token == "all")
            mechanisms.insert(make_shared<EM>(technology.getValue(), em.getValue()));
        if (token == "hci" || token == "all")
            mechanisms.insert(make_shared<HCI>(technology.getValue(), hci.getValue()));
        if (token == "tddb" || token == "all")
            mechanisms.insert(make_shared<TDDB>(technology.getValue(), tddb.getValue()));
        if (token != "all" && token != "nbti" && token != "em" && token != "hci" && token != "tddb")
            warn("ignoring unknown aging mechanism \"%s\"\n", token.c_str());
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
            cerr << "unknown unit type \"" << child.attribute("type").value()
                 << "\" for unit " << child.attribute("name").value() << endl;
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
            {
                warn("no unit failure during iteration %d\n", i);
                break;
            }

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
    cout << "Lifetime statistics for " << root->name << endl;
    cout << "Mean: " << convert_time(root->mttf(), time.getValue()) << endl;
    cout << "Standard deviation: " << convert_time(root->stdttf(), time.getValue()) << endl;
    pair<double, double> interval = root->mttf_interval(0.95);
    cout << "95\% confidence interval: [" << convert_time(interval.first, time.getValue()) << ", " << convert_time(interval.second, time.getValue()) << ']' << endl;

    if (!rates.getValue().empty())
    {
        unordered_map<string, function<double(const shared_ptr<Unit>&)>> outputs = {
            {"mttf", [&](const shared_ptr<Unit>& u){ return convert_time(u->mttf(), time.getValue()); }},
            {"failures", [](const shared_ptr<Unit>& u){ return u->ttfs.size(); }},
            {"alpha", [&](const shared_ptr<Unit>& u){ return convert_time(u->aging_rate(), time.getValue()); }}
        };
        writecsv(rates.getValue(), units, outputs);
    }
    if (!separate.getValue().empty())
    {
        unordered_map<string, function<double(const shared_ptr<Unit>&)>> outputs;
        for (const shared_ptr<FailureMechanism>& mechanism: mechanisms)
            outputs[mechanism->name] = [&](const shared_ptr<Unit>& u){ return convert_time(u->aging_rate(mechanism), time.getValue()); };
        writecsv(separate.getValue(), units, outputs);
    }
    if (!dist_dump.getValue().empty())
    {
        ofstream dist(dist_dump.getValue());
        if (dist)
        {
            dist << root->name;
            for (double ttf: root->ttfs)
                dist << ',' << convert_time(ttf, time.getValue());
            dist << endl;
            for (const shared_ptr<Unit>& unit: units)
            {
                dist << unit->name;
                for (double ttf: unit->ttfs)
                    dist << ',' << convert_time(ttf, time.getValue());
                dist << endl;
            }
        }
        else
            cerr << "error: could not write to " << dist_dump.getValue() << endl;
    }

    return 0;
}