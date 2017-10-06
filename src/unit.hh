#pragma once

#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <pugixml.hpp>
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
  protected:
    DataPoint defaults; // default values if they are missing from traces
    double peak_power;  // W
    bool fail;

    std::map<std::shared_ptr<FailureMechanism>, std::shared_ptr<ReliabilityDistribution>> reliabilities;

  public:
    Unit(const pugi::xml_node& node);
    virtual double activity(const DataPoint& data) const;
    void computeReliability(const std::vector<std::shared_ptr<FailureMechanism>>& mechanisms, std::vector<DataPoint>& trace);
    virtual double reliability(double t) const;
    virtual double mttf() const;
    virtual double mttf(const std::shared_ptr<FailureMechanism>& mechanism) const;
    bool failed() const override { return fail; }
    virtual std::ostream& dump(std::ostream& stream) const override;
};

class Group : public Component
{
  private:
    int failures;
    std::vector<std::shared_ptr<Component>> _children;

  public:
    Group(const pugi::xml_node& node, std::vector<std::shared_ptr<Unit>>& units);
    const std::vector<std::shared_ptr<Component>>& children() const { return _children; }
    double mttf() const override;
    double mttf(const std::shared_ptr<FailureMechanism>& mechanism) const override;
    bool failed() const;
    std::ostream& dump(std::ostream& ostream) const override;
};