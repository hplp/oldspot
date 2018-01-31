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
#include <set>
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

/**
 * Push a list of failed units' names in a configuration onto a stream.
 */
ostream& operator<<(ostream& os, const Unit::config_t& config)
{
    if (config.empty())
        return os << "[]";
    return os << '[' << accumulate(next(config.begin()), config.end(), *config.begin(),
                                   [](const string& a, const string& b){ return a + ',' + b; }) << ']';
}

/**
 * Get the mean of the times to failure of this Component.
 */
double Component::mttf() const
{
    if (ttfs.empty())
        return numeric_limits<double>::quiet_NaN();
    else
        return accumulate(ttfs.begin(), ttfs.end(), 0.0)/ttfs.size();
}

/**
 * Get the confidence interval on the MTTF of this Component.  Until a good implementation
 * of Student's t distribution can be found, the confidence parameter is reserved and
 * the interval is always a 95% confidence interval.
 */
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

/**
 * Push the string version of this Component onto a stream.
 */
ostream& operator<<(ostream& stream, const Component& c)
{
    return c.dump(stream);
}

// Default delimiter for parsing trace files
char Unit::delim = ',';
// config_t that specifies a "fresh" system (all units healthy)
const Unit::config_t Unit::fresh = {""};

/**
 * Check for units whose parent groups have reported failure and then mark them
 * as having failed.
 */
vector<shared_ptr<Unit>> Unit::parents_failed(const shared_ptr<Component>& root,
                                              const vector<shared_ptr<Unit>>& units)
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

/**
 * Constructor for Unit.  The new Unit reads the trace specified in the given
 * pugixml node using the given set of default values if they are missing.  Each
 * node should consist of a set of trace files that contains one trace for each
 * possible configuration in which this Unit is not failed and, optionally,
 * a redundancy specification that specifies how many redudnant copies of this Unit
 * there are and whether they are parallel or serial (i.e. if they are shadow copies
 * or take over when older ones fail).  Configurations are specified as sets of
 * names of units that have failed.
 * 
 * The ID of each unit should be unique.
 */
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
    if (traces.count(fresh) == 0)
        traces[fresh] = {{1, 1, defaults}};
    for (auto& trace: traces)
        for (DataPoint& data: trace.second)
            data.data["frequency"] *= 1e6; // Expecting MHz; convert to Hz
}

/**
 * Units don't have children, so return an empty vector.
 */
vector<shared_ptr<Component>>& Unit::children()
{
    static vector<shared_ptr<Component>> no_children;
    return no_children;
}

/**
 * Reset the unit's reliability and age to being fresh.
 */
void Unit::reset()
{
    age = 0;
    _current_reliability = 1;
    _failed = false;
    remaining = copies;
}

/**
 * Determine the configuration of failed units and set this Unit's reliability function
 * based on that.
 */
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
        config = fresh;
        cerr << "         using configuration " << config << endl;
    }
}

/**
 * Determine the next event this Enit experiences relative to the time of the previous
 * event.  Currently this only means the time at which this Unit will fail.
 */
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

/**
 * Update this Unit's reliability for the current simulation time.  This must
 * be shifted according to:
 * [1] Bolchini, C., Carminati, M., Gribaudo, M., and Miele, A. A lightweight and
 *     open-source framework for the lifetime estimation of multicore systems. ICCD 2014.
 */
void Unit::update_reliability(double dt)
{
    age += dt;
    if (!prev_config.empty())
        age -= inverse(prev_config, _current_reliability) - inverse(config, _current_reliability);
    _current_reliability = reliability(age);
}

/**
 * The activity for an unspecified type of Unit is specified directly by the trace file in an
 * "activity" column.
 */
double Unit::activity(const DataPoint& data, const shared_ptr<FailureMechanism>& mechanism) const
{
    return data.data.at("activity");
}

/**
 * Compute the reliability functions, R(t), for this Unit for all configurations.
 */
void Unit::compute_reliability(const set<shared_ptr<FailureMechanism>>& mechanisms)
{
    for (const auto& trace: traces)
    {
        for (const shared_ptr<FailureMechanism>& mechanism: mechanisms)
        {
            vector<MTTFSegment> mttfs(trace.second.size());
            for (size_t j = 0; j < trace.second.size(); j++)
            {
                double duty_cycle = min(activity(trace.second[j], mechanism), 1.0);
                double dt = j > 0 ? trace.second[j].time - trace.second[j - 1].time : trace.second[j].time;
                mttfs[j] = {dt, mechanism->timeToFailure(trace.second[j], duty_cycle)};
            }
            reliabilities[trace.first][mechanism] = mechanism->distribution(mttfs);
        }
        overall_reliabilities[trace.first] = reliabilities[trace.first].begin()->second;
        for (auto it = next(reliabilities[trace.first].begin()); it != reliabilities[trace.first].end(); ++it)
            overall_reliabilities[trace.first] *= it->second;
    }
}

/**
 * Estimate the aging rate of the unit for the given configuration if it has
 * not failed in that configuration.
 */
double Unit::aging_rate(const config_t& c) const
{
    if (failed_in_trace(c))
        return 0;
    else
        return overall_reliabilities.at(c).rate();
}

/**
 * Estimate this Unit's overall aging rate for the given failure mechanism assuming
 * a fresh system.
 */
double Unit::aging_rate(const std::shared_ptr<FailureMechanism>& mechanism) const
{
    return reliabilities.at(fresh).at(mechanism).rate();
}

/**
 * Compute this Unit's reliability at time t for configuration c.
 */
double Unit::reliability(const config_t& c, double t) const
{
    return overall_reliabilities.at(c)(t);
}

/**
 * Compute the amount of time it takes for this unit to reach reliability r with
 * configuration c.
 */
double Unit::inverse(const config_t& c, double r) const
{
    return overall_reliabilities.at(c).inverse(r);
}

/**
 * Check if this Unit is failed in the given configuration.
 */
bool Unit::failed_in_trace(const config_t& c) const
{
    return c.count(name) > 0;
}

/**
 * Set this Unit as having failed.  If there is redundancy, the amount of available
 * redudant units is decremented instead, and this Unit only fails if there are none
 * left.  If the type of redundancy is serial, then this Unit's age and reliability
 * are reset.
 */
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

/**
 * When pushed to a stream, push a Unit's name only.
 */
ostream& Unit::dump(ostream& stream) const
{
    return stream << name;
}

/**
 * Due to the complexity of a core, the activity of a Core is estimated to be the
 * fraction it is consuming of the maximum amount of power it can consume, or
 * power/peak_power.
 */
double Core::activity(const DataPoint& data, const shared_ptr<FailureMechanism>&) const
{
    return data.data.at("power")/data.data.at("peak_power");
}

/**
 * Logic activity is taken to be the number of times it is activated in a period of time
 * divided by the number of cycles that have taken place.
 * 
 * Because not all transistors in the unit are necessarily experiencing NBTI at the same
 * time, the integral of expected activity factors over all of the transistors is used instead.
 * See:
 * [2] F. Oboril and M. B. Tahoori, "ExtraTime: Modeling and analysis of
 *     wearout due to transistor aging at microarchitecture-level," in
 *     IEEE/IFIP International Conference on Dependable Systems and Networks
 *     (DSN 2012), 2012, pp. 1–12.
 */
double Logic::activity(const DataPoint& data, const shared_ptr<FailureMechanism>& mechanism) const
{
    double duty_cycle = min(data.data.at("activity")/(data.duration*data.data.at("frequency")), 1.0);
    if (mechanism->name == "NBTI")
        return 1 - duty_cycle*duty_cycle/2;
    else
        return duty_cycle;
}

/**
 * Unlike other types of units, the activity of memory units is data-dependent rather than
 * usage-dependent.  We assume here that high-order bits tend to be zero, which means their
 * degradation (particularly for NBTI) will dominate that of lower-order bits.
 */
double Memory::activity(const DataPoint& data, const shared_ptr<FailureMechanism>& mechanism) const
{
    if (mechanism->name == "HCI")
        return 0;
    else
        return 1;
}

/**
 * Constructor for a group of components.  Has a set of children that can either be other Groups
 * or Units, and is considered to be failed if enough of its children have failed.
 */
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

/**
 * A Group has failed if the number of children who have failed has passed the threshold
 * of tolerable failures.
 */
bool Group::failed() const
{
    unsigned int f = 0;
    for (const shared_ptr<Component>& child: _children)
        if (child->failed() && ++f > failures)
            return true;
    return false;
}

/**
 * Push the string representation of this Group onto a stream, which is its name followed
 * by a list of its children's names and how many failures it can tolerate.
 */
ostream& Group::dump(ostream& stream) const
{
    return stream << name << '(' << _children.size() << " children,failures=" << failures << ')';
}

} // namespace oldspot