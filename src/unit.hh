#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <pugixml.hpp>
#include <string>
#include <vector>

#include "failure.hh"
#include "reliability.hh"
#include "trace.hh"

class Component
{
  public:
    const std::string name;

    Component(const std::string _n) : name(_n) {}
    virtual double mttf() const = 0;
    virtual double mttf(const std::shared_ptr<FailureMechanism>& mechanism) const = 0;
    virtual std::ostream& dump(std::ostream& stream) const = 0;
    virtual bool failed() const = 0;

    friend std::ostream& operator<<(std::ostream& stream, const Component& c);
};

class Unit : public Component
{
  private:
    static int index;
    static std::map<uint64_t, int> trace_indices;

  protected:
    double peak_power;  // W
    bool _failed;

    std::vector<std::vector<DataPoint>> traces;
    std::vector<std::map<std::shared_ptr<FailureMechanism>, WeibullDistribution>> reliabilities;
    std::vector<WeibullDistribution> overall_reliabilities;

  public:
    static char delim;

    double current_reliability;

    static void add_configuration(uint64_t config);
    static void set_configuration(const std::vector<std::shared_ptr<Unit>>& units);

    Unit(const pugi::xml_node& node);
    void reset();
    virtual double activity(const DataPoint& data) const;
    void computeReliability(const std::vector<std::shared_ptr<FailureMechanism>>& mechanisms);
    double reliability(double t) const;
    virtual double reliability(int i, double t) const;
    double inverse(double r) const;
    virtual double inverse(int i, double r) const;
    double mttf() const override;
    virtual double mttf(int i) const;
    double mttf(const std::shared_ptr<FailureMechanism>& mechanism) const override;
    virtual double mttf(int i, const std::shared_ptr<FailureMechanism>& mechanism) const;
    bool failed() const { return _failed; }
    void failed(bool f) { _failed = f; }
    virtual std::ostream& dump(std::ostream& stream) const override;
};

class Group : public Component
{
  private:
    unsigned int failures;
    std::vector<std::shared_ptr<Component>> _children;

  public:
    Group(const pugi::xml_node& node, std::vector<std::shared_ptr<Unit>>& units);
    const std::vector<std::shared_ptr<Component>>& children() const { return _children; }
    double mttf() const override;
    double mttf(const std::shared_ptr<FailureMechanism>& mechanism) const override;
    bool failed() const;
    std::ostream& dump(std::ostream& ostream) const override;
};