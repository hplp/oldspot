#include "unit.hh"

#include <algorithm>
#include <bitset>
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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "failure.hh"
#include "reliability.hh"
#include "trace.hh"
#include "util.hh"

namespace oldspot
{

using namespace pugi;
using namespace std;

ostream& operator<<(ostream& os, const Unit::config_t& config)
{
    if (config.empty())
        return os << "[]";
    return os << '[' << accumulate(next(config.begin()), config.end(), *config.begin(),
                                   [](const string& a, const string& b){ return a + ',' + b; }) << ']';
}

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

char Unit::delim = ',';

vector<shared_ptr<Unit>> Unit::parents_failed(const shared_ptr<Component>& root, const vector<shared_ptr<Unit>>& units)
{
    vector<shared_ptr<Unit>> failed(units.begin(), units.end());
    Component::conditional_walk(root, [&](const shared_ptr<Component>& c){
        if (c->failed())
            return false;
        shared_ptr<Unit> u = dynamic_pointer_cast<Unit>(c);
        if (u)
            failed.erase(remove(failed.begin(), failed.end(), u));
        return true;
    });
    for (shared_ptr<Unit>& unit: failed)
        unit->_failed = true;
    return failed;
}

Unit::Unit(const xml_node& node, unsigned int i, unordered_map<string, double> defaults)
    : Component(node.attribute("name").value()),
      age(0), copies(1), _current_reliability(1), _failed(false), remaining(1), serial(true), config({}), prev_config({}), id(i)
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
            vector<DataPoint> trace = parseTrace(child.attribute("file").value(), delim);
            vector<string> failed_vector = split(child.attribute("failed").value(), ',');
            config_t failed(failed_vector.begin(), failed_vector.end());

            for (const auto& def: defaults)
                for (DataPoint& data: trace)
                    if (data.data.count(def.first) == 0)
                        data.data[def.first] = def.second;
            traces[failed] = trace;
        }
    }
    if (traces.count({""}) == 0)
        traces[{""}] = {{1, 1, defaults}};
    for (auto& trace: traces)
        for (DataPoint& data: trace.second)
            data.data["frequency"] *= 1e6; // Expecting MHz; convert to Hz
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

void Unit::set_configuration(const shared_ptr<Component>& root)
{
    if (_failed)
        cerr << "warning: setting configuration for failed unit " << name << endl;
    if (root->failed())
        cerr << "warning: setting configuration for failed system" << endl;

    prev_config = config;
    config.clear();
    conditional_walk(root, [&](const shared_ptr<Component>& c){
        if (c->failed())
        {
            config.insert(c->name);
            return false;
        }
        return true;
    });
    if (config.empty())
        config.insert("");

    if (traces.count(config) == 0)
    {
        cerr << "warning: can't find configuration " << config << " for " + name << endl;
        config = {""};
        cerr << "         using configuration " << config << endl;
    }
}

double Unit::get_next_event() const
{
    static random_device dev;
    static mt19937 gen(dev());
    uniform_real_distribution<double> r(0, _current_reliability);
    double next = inverse(r(gen));
    if (isinf(next))
        return numeric_limits<double>::infinity();
    return next - inverse(_current_reliability);
}

void Unit::update_reliability(double dt)
{
    age += dt;
    if (!prev_config.empty())
        age -= inverse(prev_config, _current_reliability) - inverse(config, _current_reliability);
    _current_reliability = reliability(age);
}

double Unit::activity(const DataPoint& data, const shared_ptr<FailureMechanism>& mechanism) const
{
    return data.data.at("activity");
}

void Unit::compute_reliability(const vector<shared_ptr<FailureMechanism>>& mechanisms)
{
    for (auto& trace: traces)
    {
        for (const shared_ptr<FailureMechanism> mechanism: mechanisms)
        {
            vector<MTTFSegment> mttfs(trace.second.size());
            for (size_t j = 0; j < trace.second.size(); j++)
            {
                trace.second[j].data["activity"] = min(activity(trace.second[j], mechanism), 1.0);
                double dt = j > 0 ? trace.second[j].time - trace.second[j - 1].time : trace.second[j].time;
                mttfs[j] = {dt, mechanism->timeToFailure(trace.second[j])};
            }
            reliabilities[trace.first][mechanism] = mechanism->distribution(mttfs);
        }
        overall_reliabilities[trace.first] = reliabilities[trace.first].begin()->second;
        for (auto it = next(reliabilities[trace.first].begin()); it != reliabilities[trace.first].end(); ++it)
            overall_reliabilities[trace.first] *= it->second;
    }
}

double Unit::aging_rate(const config_t& c) const
{
    if (failed_in_trace(c))
        return 0;
    else
        return overall_reliabilities.at(c).rate();
}

double Unit::reliability(const config_t& c, double t) const
{
    return overall_reliabilities.at(c)(t);
}

double Unit::inverse(const config_t& c, double r) const
{
    return overall_reliabilities.at(c).inverse(r);
}

bool Unit::failed_in_trace(const config_t& c) const
{
    return c.count(name) > 0;
}

void Unit::failure()
{
    _failed = --remaining == 0;
    if (serial)
    {
        _current_reliability = 1;
        age = 0;
        prev_config = {};
    }
}

ostream& Unit::dump(ostream& stream) const
{
    return stream << name;
}

double Core::activity(const DataPoint& data, const shared_ptr<FailureMechanism>&) const
{
    return data.data.at("power")/data.data.at("peak_power");
}

double Logic::activity(const DataPoint& data, const shared_ptr<FailureMechanism>& mechanism) const
{
    double duty_cycle = min(data.data.at("activity")/(data.duration*data.data.at("frequency")), 1.0);
    if (mechanism == NBTI::model())
        return 1 - duty_cycle*duty_cycle/2;
    else
        return duty_cycle;
}

double Memory::activity(const DataPoint& data, const shared_ptr<FailureMechanism>& mechanism) const
{
    if (mechanism == HCI::model())
        return 0;
    else
        return 1;
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
    unsigned int f = 0;
    for (const shared_ptr<Component>& child: _children)
        if (child->failed() && ++f > failures)
            return true;
    return false;
}

ostream& Group::dump(ostream& stream) const
{
    return stream << name << '(' << _children.size() << " children,failures=" << failures << ')';
}

} // namespace oldspot