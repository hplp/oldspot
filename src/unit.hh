#pragma once

#include <bitset>
#include <cstdint>
#include <limits>
#include <memory>
#include <ostream>
#include <pugixml.hpp>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "failure.hh"
#include "reliability.hh"
#include "trace.hh"

namespace std
{

/**
 * Hash functor that allows unordered_set<string> to be used as a map key.
 */
template <>
class hash<unordered_set<string>>
{
  public:
    size_t
    operator()(const unordered_set<string>& strs) const
    {
        size_t result = 58271;
        for (const string& s: strs)
            result ^= strs.hash_function()(s);
        return result;
    }
};

}

namespace oldspot
{

/**
 * Component in the system to be simulated.  This can either be a Group, which
 * contains other components and has its failure state depend on the failure
 * states of its children, or a Unit, which has no children and has its failure
 * state depend on reliability calculation.
 */
class Component
{
  public:
    /**
     * Perform a function the given component and each of its children in a
     * prefix depth-first traversal.
     */
    template<typename function> static void
    walk(const std::shared_ptr<Component>& root, function&& op)
    {
        using namespace std;
    
        stack<shared_ptr<Component>> components;
        components.push(root);
        while (!components.empty())
        {
            shared_ptr<Component> component = components.top();
            components.pop();
            op(component);
            for (const shared_ptr<Component>& child: component->children())
                components.push(child);
        }
    }

    /**
     * Perform a function the given component and each of its children in a
     * prefix depth-first traversal.  A component's children are only traversed
     * if the function performed on that component returns a value that evaluates
     * to true.
     */
    template<typename function> static void
    conditional_walk(const std::shared_ptr<Component>& root, function&& op)
    {
        using namespace std;
    
        stack<shared_ptr<Component>> components;
        components.push(root);
        while (!components.empty())
        {
            shared_ptr<Component> component = components.top();
            components.pop();
            if (op(component))
                for (const shared_ptr<Component>& child: component->children())
                    components.push(child);
        }
    }

    const std::string name;
    std::vector<double> ttfs;

    Component(const std::string _n) : name(_n) {}
    virtual std::vector<std::shared_ptr<Component>>& children() = 0;
    virtual double mttf() const;
    virtual std::pair<double, double> mttf_interval(double confidence=0.95) const;
    virtual double aging_rate() const { return std::numeric_limits<double>::quiet_NaN(); }
    virtual bool failed() const = 0;

    virtual std::ostream& dump(std::ostream& stream) const = 0;
    friend std::ostream& operator<<(std::ostream& stream, const Component& c);
};

/**
 * A unit in the system, represented as a leaf node in the failure dependency
 * graph.  Each unit is associated with a trace of power, performance, temperature,
 * etc. that affects the rate at which its reliability degrades.  Each unit requires
 * one of these traces for each healthy configuration of the system except for ones
 * on which the unit has failed.
 */
class Unit : public Component
{
  public:
    typedef std::unordered_set<std::string> config_t;
    typedef enum { HEALTHY, ASLEEP, FAILED } state_t;

  private:
    double _area;
    double age;
    int copies;
    double _current_reliability;
    state_t _state;
    int remaining;
    bool serial;
    config_t config;
    config_t prev_config;
    double updowntime;

  protected:
    std::unordered_map<config_t, std::vector<DataPoint>> traces;
    std::unordered_map<config_t, std::pair<double, double>> power_gating; // uptime, downtime
    std::unordered_map<config_t, std::unordered_map<std::shared_ptr<FailureMechanism>, WeibullDistribution>> reliabilities;
    std::unordered_map<config_t, WeibullDistribution> overall_reliabilities;

  public:
    static char delim;
    static const config_t fresh;

    static std::vector<std::shared_ptr<Unit>> parents_failed(const std::shared_ptr<Component>& root, const std::vector<std::shared_ptr<Unit>>& units);

    double area() const { return _area; }
    const unsigned int id;

    Unit(const pugi::xml_node& node, unsigned int i, std::unordered_map<std::string, double> defaults={});
    std::vector<std::shared_ptr<Component>>& children() override;
    void reset();
    void set_configuration(const std::shared_ptr<Component>& root);

    std::pair<double, state_t> get_next_event() const;
    void update_reliability(double dt);
    double current_reliability() const { return _current_reliability; }

    virtual double activity(const DataPoint& data, const std::shared_ptr<FailureMechanism>& mechanism) const;
    void compute_reliability(const std::set<std::shared_ptr<FailureMechanism>>& mechanisms);

    double aging_rate(const config_t& c) const;
    double aging_rate() const override { return aging_rate(fresh); }
    double aging_rate(const std::shared_ptr<FailureMechanism>& mechanism) const;

    virtual double reliability(const config_t& c, double t) const;
    double reliability(double t) const { return reliability(config, t); }

    virtual double inverse(const config_t& c, double r) const;
    double inverse(double r) const { return inverse(config, r); }

    bool failed_in_trace(const config_t& c) const;
    bool failed() const { return _state == FAILED; }
    void failure();

    virtual std::ostream& dump(std::ostream& stream) const override;
};

std::ostream& operator<<(std::ostream& os, const Unit::config_t& config);

/**
 * Unit that represents an entire core.  The average activity factor of a core
 * is estimated as its current power consumption divided by its peak power
 * consumption, or the maximum amount of power it can draw.
 */
class Core : public Unit
{
  public:
    Core(const pugi::xml_node& node, unsigned int i)
        : Unit(node, i, {{"power", 1}, {"peak_power", 1}}) {}
    double activity(const DataPoint& data, const std::shared_ptr<FailureMechanism>&) const override;
};

/**
 * Unit that represents a block consisting primarily of logic gates.
 */
class Logic : public Unit
{
  public:
    Logic(const pugi::xml_node& node, unsigned int i) : Unit(node, i) {}
    double activity(const DataPoint& data, const std::shared_ptr<FailureMechanism>&) const override;
};

/**
 * Unit that represents a block consisting primarily of memory elements.  Unlike
 * other types of units, its activity is data-dependent rather than usage-
 * dependent.
 */
class Memory : public Unit
{
  public:
    Memory(const pugi::xml_node& node, unsigned int i) : Unit(node, i) {}
    double activity(const DataPoint& data, const std::shared_ptr<FailureMechanism>& mechanism) const override;
};

/**
 * Group component that has children which are either other groups or Units.
 * A group has a threshold that defines how many failures in its children it
 * can tolerate before it itself has failed.
 */
class Group : public Component
{
  private:
    unsigned int failures;
    std::vector<std::shared_ptr<Component>> _children;

  public:
    Group(const pugi::xml_node& node, std::vector<std::shared_ptr<Unit>>& units);
    std::vector<std::shared_ptr<Component>>& children() override { return _children; }
    bool failed() const;
    std::ostream& dump(std::ostream& ostream) const override;
};

} // namespace oldspot