#pragma once

#include <bitset>
#include <cstdint>
#include <limits>
#include <memory>
#include <ostream>
#include <pugixml.hpp>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "failure.hh"
#include "reliability.hh"
#include "trace.hh"

namespace oldspot
{

class Component
{
  public:
    template<typename function>
    static void walk(const std::shared_ptr<Component>& root, function&& op)
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

    const std::string name;
    std::vector<double> ttfs;

    Component(const std::string _n) : name(_n) {}
    virtual std::vector<std::shared_ptr<Component>>& children() = 0;
    virtual double mttf() const;
    virtual std::pair<double, double> mttf_interval(double confidence=0.95) const;
    virtual bool failed() const = 0;

    virtual std::ostream& dump(std::ostream& stream) const = 0;
    friend std::ostream& operator<<(std::ostream& stream, const Component& c);
};

class Unit : public Component
{
  private:
    static std::unordered_map<std::bitset<64>, int> trace_indices;

    double age;
    int copies;
    double _current_reliability;
    bool _failed;
    int remaining;
    bool serial;
    int index;
    int prev_index;

  protected:
    std::vector<std::vector<DataPoint>> traces;
    std::vector<std::unordered_map<std::shared_ptr<FailureMechanism>, WeibullDistribution>> reliabilities;
    std::vector<WeibullDistribution> overall_reliabilities;

  public:
    static char delim;

    static std::vector<std::vector<std::shared_ptr<Unit>>> init_configurations(const std::shared_ptr<Component>& root, std::vector<std::shared_ptr<Unit>>& units);
    static void set_configuration(const std::vector<std::shared_ptr<Unit>>& units);
    static size_t configurations() { return trace_indices.size(); }

    const unsigned int id;

    Unit(const pugi::xml_node& node, unsigned int i, std::unordered_map<std::string, double> defaults={});
    std::vector<std::shared_ptr<Component>>& children() override;
    void reset();
    double get_next_event() const;
    void update_reliability(double dt);
    double current_reliability() const { return _current_reliability; }

    virtual double activity(const DataPoint& data, const std::shared_ptr<FailureMechanism>& mechanism) const;
    void computeReliability(const std::vector<std::shared_ptr<FailureMechanism>>& mechanisms);

    double aging_rate(int i) const;
    double aging_rate() const { return aging_rate(index); }

    double reliability(double t) const { return reliability(index, t); }
    virtual double reliability(int i, double t) const;

    double inverse(double r) const { return inverse(index, r); }
    virtual double inverse(int i, double r) const;

    bool failed_in_trace(int i) const;
    bool failed() const { return _failed; }
    void failure();

    virtual std::ostream& dump(std::ostream& stream) const override;
};

class Core : public Unit
{
  public:
    Core(const pugi::xml_node& node, unsigned int i)
        : Unit(node, i, {{"power", 1}, {"peak_power", 1}}) {}
    double activity(const DataPoint& data, const std::shared_ptr<FailureMechanism>&) const override;
};

class Logic : public Unit
{
  public:
    Logic(const pugi::xml_node& node, unsigned int i) : Unit(node, i) {}
    double activity(const DataPoint& data, const std::shared_ptr<FailureMechanism>&) const override;
};

class Memory : public Unit
{
  public:
    Memory(const pugi::xml_node& node, unsigned int i) : Unit(node, i) {}
    // This is data-dependent rather than usage-dependent, but we assume that
    // high-order bits tend to be zero, which gives an activity factor of 1
    // in SRAM (we also assume that the aging of the SRAM dominates that of
    // the addressing logic)
    double activity(const DataPoint& data, const std::shared_ptr<FailureMechanism>& mechanism) const override;
};

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