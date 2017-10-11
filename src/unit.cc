#include "unit.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <pugixml.hpp>
#include <utility>
#include <vector>

#include "failure.hh"
#include "interp.hh"
#include "reliability.hh"
#include "trace.hh"

using namespace pugi;
using namespace std;

ostream& operator<<(ostream& stream, const Component& c)
{
    return c.dump(stream);
}

int Unit::index = 0;

map<uint64_t, int> Unit::trace_indices;

char Unit::delim = ',';

void Unit::add_configuration(uint64_t config)
{
    trace_indices[config] = index++;
}

void Unit::set_configuration(const vector<shared_ptr<Unit>>& units)
{
    uint64_t config = accumulate(units.begin(), units.end(), 0,
                                 [](uint64_t a, const shared_ptr<Unit>& b){ return (a << 1) | (b->failed() ? 0 : 1); });
    index = trace_indices.at(config);
}

Unit::Unit(const xml_node& node)
    : Component(node.attribute("name").value()), peak_power(0), _failed(false), current_reliability(1)
{
    map<string, double> defaults = {{"vdd", 1}, {"temperature", 350}, {"frequency", 1000}};
    for (const xml_node& def: node.children("default"))
    {
        if (def.attribute("vdd"))
            defaults["vdd"] = def.attribute("vdd").as_double();
        if (def.attribute("temperature"))
            defaults["temperature"] = def.attribute("temperature").as_double();
        if (def.attribute("frequency"))
            defaults["frequency"] = def.attribute("frequency").as_double();
        if (def.attribute("peak_power"))
            peak_power = def.attribute("peak_power").as_double();
    }

    if (node.child("trace"))
    {
        for (const xml_node& child: node.children("trace"))
        {
            traces.push_back(parseTrace(child.attribute("file").value(), delim));
            for (const auto& def: defaults)
                for (DataPoint& data: traces.back())
                    if (data.data.count(def.first) == 0)
                        data.data[def.first] = def.second;
        }
    }
    else
        traces.push_back({{1, defaults}});
    for (vector<DataPoint>& trace: traces)
        for (DataPoint& data: trace)
            data.data["frequency"] *= 1e6; // Expecting MHz; convert to Hz

    reliabilities.resize(traces.size());
    overall_reliabilities.resize(traces.size());
}

void Unit::reset()
{
    _failed = false;
    current_reliability = 1;
}

double Unit::activity(const DataPoint& data) const
{
    // (cycles where unit is active)/(cycles of time step)
    // (runtime power)/(peak power)
    // For NBTI in an SRAM, this is data-dependent rather than usage-dependent
    return data.data.at("activity");
}

void Unit::computeReliability(const vector<shared_ptr<FailureMechanism>>& mechanisms)
{
    for (size_t i = 0; i < traces.size(); i++)
    {
        for (const shared_ptr<FailureMechanism> mechanism: mechanisms)
        {
            vector<MTTFSegment> mttfs(traces[i].size());
            for (size_t j = 0; j < traces[i].size(); j++)
            {
                traces[i][j].data["activity"] = activity(traces[i][j]);
                double dt = j > 0 ? traces[i][j].time - traces[i][j - 1].time : traces[i][j].time;
                mttfs[j] = {dt, mechanism->timeToFailure(traces[i][j])};
            }
            reliabilities[i][mechanism] = mechanism->distribution(mttfs);
        }
        overall_reliabilities[i] = reliabilities[i].begin()->second;
        for (auto it = next(reliabilities[i].begin()); it != reliabilities[i].end(); ++it)
            overall_reliabilities[i] *= it->second;
    }
}

double Unit::reliability(double t) const
{
    return reliability(index, t);
}

double Unit::reliability(int i, double t) const
{
    return overall_reliabilities[i].reliability(t);
}

double Unit::inverse(double r) const
{
    return inverse(index, r);
}

double Unit::inverse(int i, double r) const
{
    return overall_reliabilities[i].inverse(r);
}

double Unit::mttf() const
{
    return mttf(index);
}

double Unit::mttf(int i) const
{
    return overall_reliabilities[i].mttf();
}

double Unit::mttf(const shared_ptr<FailureMechanism>& mechanism) const
{
    return mttf(index, mechanism);
}

double Unit::mttf(int i, const shared_ptr<FailureMechanism>& mechanism) const
{
    return reliabilities[i].at(mechanism).mttf();
}

ostream& Unit::dump(ostream& stream) const
{
    return stream << name;
}

Group::Group(const xml_node& node, vector<shared_ptr<Unit>>& units)
    : Component(node.attribute("name").value()), failures(node.attribute("failures").as_int())
{
    for (const xml_node& child: node.children())
    {
        if (strcmp(child.name(), "group") == 0)
            _children.push_back(make_shared<Group>(child, units));
        else if (strcmp(child.name(), "unit") == 0)
        {
            shared_ptr<Unit> unit = make_shared<Unit>(child);
            units.push_back(unit);
            _children.push_back(unit);
        }
        else
        {
            cerr << "unknown component type " << child.name() << endl;
            exit(1);
        }
    }
}

double Group::mttf() const
{
    return 0;
}

double Group::mttf(const shared_ptr<FailureMechanism>& mechanism) const
{
    return 0;
}

bool Group::failed() const
{
    unsigned int f = count_if(_children.begin(), _children.end(),
                              [](const shared_ptr<Component>& c){ return c->failed(); });
    return f > failures;
}

ostream& Group::dump(ostream& stream) const
{
    return stream << name << '(' << _children.size() << " children,failures=" << failures << ')';
}