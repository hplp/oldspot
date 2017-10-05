#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <pugixml.hpp>
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

    try
    {
        CmdLine cmd("Compute the reliability distribution of a chip", ' ', "0.1");
        ValueArg<char> delimiter("", "trace-delimiter", "One-character delimiter for data in input trace files (default: ,)", false, ',', "delim", cmd);
        ValueArg<string> ptrace("p", "power", "File containing power traces for all units (W)", false, "", "filename", cmd);
        ValueArg<string> ftrace("f", "frequency", "File containing frequency traces for all units (MHz)", false, "", "filename", cmd);
        ValueArg<string> ttrace("T", "temperature", "File containing temperature traces for all units (K)", false, "", "filename", cmd);
        ValueArg<string> vtrace("v", "vdd", "File containing voltage traces for all units (V)", false, "", "filename", cmd);
        ValueArg<string> atrace("a", "activity", "File containing activity traces for all units", false, "", "filename", cmd);
        ValueArg<string> config("c", "chip-config", "File containing chip configuration", true, "", "filename", cmd);
        cmd.parse(argc, argv);

        xml_document doc;
        xml_parse_result result = doc.load_file(config.getValue().c_str());
        if (!result)
        {
            cerr << config.getValue() << ": " << result.description()
                 << " at offset " << result.offset << endl;
            return 1;
        }

        vector<shared_ptr<Unit>> units;
        Group root(doc, units);

        trace_t activity, voltage, temperature, frequency, power;
        if (!atrace.getValue().empty())
            activity = parseTrace(atrace.getValue(), delimiter.getValue());
        if (!vtrace.getValue().empty())
            voltage = parseTrace(vtrace.getValue(), delimiter.getValue());
        if (!ttrace.getValue().empty())
            temperature = parseTrace(ttrace.getValue(), delimiter.getValue());
        if (!ftrace.getValue().empty())
            frequency = parseTrace(ftrace.getValue(), delimiter.getValue());
        if (!ptrace.getValue().empty())
            power = parseTrace(ptrace.getValue(), delimiter.getValue());
        unordered_set<string> names(units.size());
        for (const shared_ptr<Unit>& unit: units)
            names.insert(unit->name);
        map<string, vector<DataPoint>> traces = collectTraces(names, activity, voltage, temperature, frequency, power);

        vector<shared_ptr<FailureMechanism>> mechanisms = {NBTI::model()};
        for (const shared_ptr<Unit>& unit: units)
            unit->computeReliability(mechanisms, traces[unit->name]);
        // Compute reliabilities
        // Monte Carlo sim to get overall failure distribution
    }
    catch (ArgException& e)
    {
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
        return 1;
    }
    return 0;
}