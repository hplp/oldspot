#pragma once

#include <cstdlib>
#include <deque>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

struct DataPoint
{
    static size_t size() { return 5; }

    double duration;    // s
    double activity;
    double vdd;         // V
    double temperature; // K
    double frequency;   // Hz
    double power;       // W

    DataPoint(double _d,
              double _a=std::numeric_limits<double>::signaling_NaN(),
              double _v=std::numeric_limits<double>::signaling_NaN(),
              double _t=std::numeric_limits<double>::signaling_NaN(),
              double _f=std::numeric_limits<double>::signaling_NaN(),
              double _p=std::numeric_limits<double>::signaling_NaN())
        : duration(_d), activity(_a), vdd(_v), temperature(_t), frequency(_f), power(_p)
    {}

    double& operator[](int i);

    friend std::ostream& operator<<(std::ostream& stream, const DataPoint& point);
};

typedef std::map<std::string, std::deque<std::pair<double, double>>> trace_t;

trace_t parseTrace(const std::string fname, char delimiter=',');
std::map<std::string, std::vector<DataPoint>> collectTraces(const std::unordered_set<std::string>& units,
                                                            trace_t& atrace,
                                                            trace_t& vtrace,
                                                            trace_t& ttrace,
                                                            trace_t& ftrace,
                                                            trace_t& ptrace);