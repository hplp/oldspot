#include "unit.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <pugixml.hpp>
#include <utility>
#include <vector>

#include "failure.hh"
#include "reliability.hh"
#include "trace.hh"

using namespace pugi;
using namespace std;

ostream& operator<<(ostream& stream, const Component& c)
{
    return c.dump(stream);
}

Unit::Unit(const xml_node& node)
    : Component(node.attribute("name").value()), defaults({0, 0, 1, 350, 1e9, 0}), peak_power(0)
{
    for (const xml_node& def: node.children("default"))
    {
        if (def.attribute("vdd"))
            defaults.vdd = def.attribute("vdd").as_double();
        if (def.attribute("temperature"))
            defaults.temperature = def.attribute("temperature").as_double();
        if (def.attribute("frequency"))
            defaults.frequency = def.attribute("frequency").as_double();
        if (def.attribute("peak_power"))
            peak_power = def.attribute("peak_power").as_double();
    }
}

double Unit::activity(const DataPoint& data) const
{
    // (cycles where unit is active)/(cycles of time step)
    // (runtime power)/(peak power)
    // For NBTI in an SRAM, this is data-dependent rather than usage-dependent
    return data.activity;
}

void Unit::computeReliability(const vector<shared_ptr<FailureMechanism>>& mechanisms, vector<DataPoint>& trace)
{
    for (size_t j = 0; j < DataPoint::size(); j++)
    {
        for (size_t i = 1; i < trace.size(); i++)
            if (isnan(trace[i][j]))
                trace[i][j] = trace[i - 1][j];
        if (isnan(trace.back()[j]))
            trace.back()[j] = defaults[j];
        for (int i = trace.size() - 2; i >= 0; i--)
            if (isnan(trace[i][j]))
                trace[i][j] = trace[i + 1][j];
    }

    for (const shared_ptr<FailureMechanism> mechanism: mechanisms)
    {
        vector<MTTFSegment> mttfs(trace.size());
        for (size_t i = 0; i < trace.size(); i++)
        {
            trace[i].activity = activity(trace[i]);
            mttfs[i] = {trace[i].duration, mechanism->timeToFailure(trace[i])};
        }
        reliabilities[mechanism] = mechanism->distribution(mttfs);
    }

    cout << name << endl;
    for (const DataPoint& point: trace)
        cout << '\t' << point.duration << ": " << point.activity << '\t' << point.vdd << '\t' << point.temperature << '\t' << point.frequency << '\t' << point.power << endl;
    cout << "\tMTTF:" << endl;
    for (const auto& mechanism: mechanisms)
        cout << "\t\t" << mechanism->name << ": " << mttf(mechanism)/(3600*24*365) << endl;
    cout << "\t\tOverall: " << mttf()/(3600*24*365) << endl;
    cout << endl;
}

double Unit::reliability(double t) const
{
    double reliability = 1;
    for (const auto& r: reliabilities)
        reliability *= r.second->reliability(t);
    return reliability;
}

double Unit::mttf() const
{
    double t = numeric_limits<double>::max();
    for (const auto& r: reliabilities)
        t = min(t, r.second->mttf());
    return t;
}

double Unit::mttf(const shared_ptr<FailureMechanism>& mechanism) const
{
    return reliabilities.at(mechanism)->mttf();
}

ostream& Unit::dump(ostream& stream) const
{
    return stream << name << "(vdd=" << defaults.vdd
                          << ",T=" << defaults.temperature
                          << ",f=" << defaults.frequency
                          << ",P=" << peak_power << ')';
}

Group::Group(const xml_node& node, vector<shared_ptr<Unit>>& units)
    : Component(node.attribute("name").value()), failures(0)
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

ostream& Group::dump(ostream& stream) const
{
    return stream << name << '(' << _children.size() << " children,failures=" << failures << ')';
}