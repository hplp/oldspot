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
#include <random>
#include <utility>
#include <vector>

#include "failure.hh"
#include "reliability.hh"
#include "trace.hh"

namespace oldspot
{

using namespace pugi;
using namespace std;

double Component::mttf() const
{
    if (ttfs.empty())
        return numeric_limits<double>::quiet_NaN();
    else
        return accumulate(ttfs.begin(), ttfs.end(), 0.0)/ttfs.size();
}

// confidence is reserved until a good implementation of Student's t distribution can be found
// right now, the confidence interval is always 95%
pair<double, double> Component::mttf_interval(double confidence) const
{
    if (ttfs.size() <= 1)
        return {numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()};
    else
    {
        double mean = mttf();
        double s = sqrt(accumulate(ttfs.begin(), ttfs.end(), 0.0,
                                   [&](double a, double b){ return a + pow(b - mean, 2); })/(ttfs.size() - 1));
        return {mean - 1.96*s/sqrt(ttfs.size()), mean + 1.96*s/sqrt(ttfs.size())};
    }
}

ostream& operator<<(ostream& stream, const Component& c)
{
    return c.dump(stream);
}

map<vector<bool>, int> Unit::trace_indices;
char Unit::delim = ',';

vector<bool> operator++(vector<bool>& n, int)
{
    vector<bool> r = n;
    for (size_t i = 0; i < n.size(); i++)
    {
        if (!n[i])
        {
            n[i] = true;
            return r;
        }
        n[i] = false;
    }
    return r;
}

void Unit::init_configurations(const shared_ptr<Component>& root, vector<shared_ptr<Unit>>& units)
{
    int index = 0;
    for (vector<bool> failed(units.size()); !all_of(failed.begin(), failed.end(), [](bool b){ return b; }); failed++)
    {
        for (const shared_ptr<Unit>& unit: units)
            unit->_failed = failed[unit->id];
        if (!root->failed())
            trace_indices[failed] = index++;
    }
    // A system where failed is all true must be failed
}

vector<vector<shared_ptr<Unit>>> Unit::valid_configurations(const shared_ptr<Component>& root, vector<shared_ptr<Unit>>& units)
{
    vector<vector<shared_ptr<Unit>>> names;
    for (vector<bool> failed(units.size()); !all_of(failed.begin(), failed.end(), [](bool b){ return b; }); failed++)
    {
        for (const shared_ptr<Unit>& unit: units)
            unit->_failed = failed[unit->id];
        if (!root->failed())
        {
            vector<shared_ptr<Unit>> functional;
            for (shared_ptr<Unit>& unit: units)
                if (!unit->failed())
                    functional.push_back(unit);
            names.push_back(functional);
        }
    }
    return names;
}

void Unit::set_configuration(const vector<shared_ptr<Unit>>& units)
{
    vector<bool> config(units.size());
    for (const shared_ptr<Unit>& unit: units)
        config[unit->id] = unit->failed();
    for (const shared_ptr<Unit>& unit: units)
    {
        unit->prev_index = unit->index;
        unit->index = trace_indices.at(config);
    }
}

Unit::Unit(const xml_node& node, unsigned int i, map<string, double> defaults)
    : Component(node.attribute("name").value()),
      age(0), copies(1), _current_reliability(1), _failed(false), remaining(1), serial(true), index(-1), prev_index(-1), id(i)
{
    if (defaults.count("vdd") == 0)
        defaults["vdd"] = 1;
    if (defaults.count("temperature") == 0)
        defaults["temperature"] = 350;
    if (defaults.count("frequency") == 0)
        defaults["frequency"] = 1000;
    if (defaults.count("activity") == 0)
        defaults["activity"] = 0;

    for (const xml_node& def: node.children("default"))
        for (auto& value: defaults)
            if (def.attribute(value.first.c_str()))
                value.second = def.attribute(value.first.c_str()).as_double();

    if (node.child("redundancy"))
    {
        const xml_node& redundancy = node.child("redundancy");
        serial = strcmp(redundancy.attribute("type").value(), "serial") == 0;
        copies = remaining = redundancy.attribute("count").as_int();
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
        traces.push_back({{1, 1, defaults}});
    for (vector<DataPoint>& trace: traces)
        for (DataPoint& data: trace)
            data.data["frequency"] *= 1e6; // Expecting MHz; convert to Hz

    reliabilities.resize(traces.size());
    overall_reliabilities.resize(traces.size());
}

vector<shared_ptr<Component>>& Unit::children()
{
    static vector<shared_ptr<Component>> no_children;
    return no_children;
}

void Unit::reset()
{
    age = 0;
    _current_reliability = 1;
    _failed = false;
    remaining = copies;
}

double Unit::get_next_event() const
{
    static random_device dev;
    static mt19937 gen(dev());
    uniform_real_distribution<double> r(0, _current_reliability);
    return inverse(r(gen)) - inverse(_current_reliability);
}

void Unit::update_reliability(double dt)
{
    age += dt;
    if (prev_index >= 0)
        age -= inverse(prev_index, _current_reliability) - inverse(index, _current_reliability);
    _current_reliability = reliability(age);
}

double Unit::activity(const DataPoint& data) const
{
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

double Unit::aging_rate(int i) const
{
    if (failed_in_trace(i))
        return 0;
    else
        return overall_reliabilities[i].rate();
}

double Unit::reliability(int i, double t) const
{
    return overall_reliabilities[i](t);
}

double Unit::inverse(int i, double r) const
{
    return overall_reliabilities[i].inverse(r);
}

bool Unit::failed_in_trace(int i) const
{
    auto result = find_if(trace_indices.begin(), trace_indices.end(),
                          [&](pair<vector<bool>, int> a){ return a.second == i; });
    if (result != trace_indices.end())
        return result->first[id];
    else
        return false;
}

void Unit::failure()
{
    _failed = --remaining == 0;
    if (serial)
    {
        _current_reliability = 1;
        age = 0;
        prev_index = -1;
    }
}

ostream& Unit::dump(ostream& stream) const
{
    return stream << name;
}

double Core::activity(const DataPoint& data) const
{
    return data.data.at("power")/data.data.at("peak_power");
}

double Logic::activity(const DataPoint& data) const
{
    return data.data.at("activity")/(data.duration*data.data.at("frequency"));
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
            string n = child.attribute("name").value();
            auto unit = find_if(units.begin(), units.end(),
                                [&](const shared_ptr<Unit>& u){ return u->name == n; });
            if (unit != units.end())
                _children.push_back(*unit);
        }
        else
        {
            cerr << "unknown component type " << child.name() << endl;
            exit(1);
        }
    }
}

bool Group::failed() const
{
    unsigned int f = count_if(_children.begin(), _children.end(),
                              [](const shared_ptr<Component>& c){ return c->failed(); });
    return f > failures || f >= _children.size();
}

ostream& Group::dump(ostream& stream) const
{
    return stream << name << '(' << _children.size() << " children,failures=" << failures << ')';
}

} // namespace oldspot