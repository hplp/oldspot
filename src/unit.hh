#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <pugixml.hpp>
#include <stack>
#include <string>
#include <vector>

#include "failure.hh"
#include "reliability.hh"
#include "trace.hh"

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

    Component(const std::string _n, size_t n=1) : name(_n) { ttfs.resize(n); }
    virtual std::vector<std::shared_ptr<Component>>& children() = 0;
    virtual double mttf() const;
    virtual bool failed() const = 0;

    virtual std::ostream& dump(std::ostream& stream) const = 0;
    friend std::ostream& operator<<(std::ostream& stream, const Component& c);
};

class Unit : public Component
{
  private:
    static int index;
    static std::map<uint64_t, int> trace_indices;

  protected:
    bool _failed;

    std::vector<std::vector<DataPoint>> traces;
    std::vector<std::map<std::shared_ptr<FailureMechanism>, WeibullDistribution>> reliabilities;
    std::vector<WeibullDistribution> overall_reliabilities;

  public:
    static char delim;

    static void add_configuration(uint64_t config);
    static void set_configuration(const std::vector<std::shared_ptr<Unit>>& units);
    static size_t configurations() { return trace_indices.size(); }

    const unsigned int id;
    double current_reliability;

    Unit(const pugi::xml_node& node, unsigned int i, size_t n=1,
         std::map<std::string, double> defaults={{"vdd", 1}, {"temperature", 350}, {"frequency", 1000}});
    std::vector<std::shared_ptr<Component>>& children() override;
    void reset();
    virtual double activity(const DataPoint& data) const;
    void computeReliability(const std::vector<std::shared_ptr<FailureMechanism>>& mechanisms);
    double aging_rate(int i) const;
    double aging_rate() const { return aging_rate(index); }
    double reliability(double t) const { return reliability(index, t); }
    virtual double reliability(int i, double t) const;
    double inverse(double r) const { return inverse(index, r); }
    virtual double inverse(int i, double r) const;
    bool failed_in_trace(int i) const;
    bool failed() const { return _failed; }
    void failed(bool f) { _failed = f; }
    virtual std::ostream& dump(std::ostream& stream) const override;
};

class Core : public Unit
{
  public:
    Core(const pugi::xml_node& node, unsigned int i, size_t n=1,
         std::map<std::string, double> defaults={{"vdd", 1}, {"temperature", 350}, {"frequency", 1000}, {"power", 1}, {"peak_power", 1}})
        : Unit(node, i, n, defaults) {}
    double activity(const DataPoint& data) const override;
};

class Logic : public Unit
{
  public:
    Logic(const pugi::xml_node& node, unsigned int i, size_t n=1) : Unit(node, i, n) {}
    double activity(const DataPoint& data) const override;
};

class Group : public Component
{
  private:
    unsigned int failures;
    std::vector<std::shared_ptr<Component>> _children;

  public:
    Group(const pugi::xml_node& node, std::vector<std::shared_ptr<Unit>>& units, size_t n=1);
    std::vector<std::shared_ptr<Component>>& children() override { return _children; }
    bool failed() const;
    std::ostream& dump(std::ostream& ostream) const override;
};