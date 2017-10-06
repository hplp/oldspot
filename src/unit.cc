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

char Unit::delim = ',';

Unit::Unit(const xml_node& node)
    : Component(node.attribute("name").value()), peak_power(0)
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
    for (const shared_ptr<FailureMechanism> mechanism: mechanisms)
    {
        for (size_t i = 0; i < traces.size(); i++)
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
    }

    cout << name << ": " << traces.size() << " trace(s)" << endl;
    for (size_t i = 0; i < traces.size(); i++)
    {
        cout << "\tTrace " << i << endl;
        for (const auto& mechanism: mechanisms)
            cout << "\t\t" << mechanism->name << ": " << mttf(mechanism, i)/(3600*24*365) << endl;
        cout << "\t\tOverall: " << mttf(i)/(3600*24*365) << endl;
        cout << endl;
    }
}

double Unit::reliability(double t, int i) const
{
    double reliability = 1;
    for (const auto& r: reliabilities[i])
        reliability *= r.second->reliability(t);
    return reliability;
}

double Unit::mttf() const
{
    return mttf(0);
}

double Unit::mttf(int i) const
{
    double t = numeric_limits<double>::max();
    for (const auto& r: reliabilities[i])
        t = min(t, r.second->mttf());
    return t;
}

double Unit::mttf(const shared_ptr<FailureMechanism>& mechanism) const
{
    return mttf(0);
}

double Unit::mttf(const shared_ptr<FailureMechanism>& mechanism, int i) const
{
    return reliabilities[i].at(mechanism)->mttf();
}

ostream& Unit::dump(ostream& stream) const
{
    return stream << name;
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

bool Group::failed() const
{
    return count_if(_children.begin(), _children.end(), [](const shared_ptr<Component>& c){ return c->failed(); }) > failures;
}

ostream& Group::dump(ostream& stream) const
{
    return stream << name << '(' << _children.size() << " children,failures=" << failures << ')';
}