#include "trace.hh"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

double& DataPoint::operator[](int i)
{
    switch (i)
    {
      case 0:
        return activity;
      case 1:
        return vdd;
      case 2:
        return temperature;
      case 3:
        return frequency;
      case 4:
        return power;
      default:
        throw out_of_range("DataPoint index out of range: " + to_string(i));
    }
}

ostream& operator<<(ostream& stream, const DataPoint& point)
{
    return stream << point.duration << ":(" << point.activity << ','
                                            << point.vdd << ','
                                            << point.temperature << ','
                                            << point.frequency << ','
                                            << point.power << ')';
}

trace_t parseTrace(const string fname, char delimiter)
{
    ifstream file(fname);
    if (!file)
    {
        cerr << fname << ": unable to open file" << endl;
        exit(1);
    }

    string line, token;
    istringstream stream;
    vector<string> units;
    trace_t trace;
    
    // Parse unit names
    getline(file, line); // Get line containing unit names
    stream.str(line);
    stream.clear();
    getline(stream, token, delimiter); // Ignore header of time column
    while (getline(stream, token, delimiter))
        units.push_back(token);

    // Parse times (first column) and values
    while (getline(file, line))
    {
        stream.str(line);
        stream.clear();
        getline(stream, token, delimiter);
        double time = stod(token); // First column should be time
        for (const string& unit: units)
        {
            getline(stream, token, delimiter);
            trace[unit].push_back({time, stod(token)});
        }
    }

    return trace;
}

map<string, vector<DataPoint>> collectTraces(const unordered_set<string>& units,
                                             trace_t& atrace,
                                             trace_t& vtrace,
                                             trace_t& ttrace,
                                             trace_t& ftrace,
                                             trace_t& ptrace)
{
    pair<double, double> missing = {numeric_limits<double>::max(), numeric_limits<double>::signaling_NaN()};
    map<string, vector<DataPoint>> trace;
    for (const string& unit: units)
    {
        double prev_time = 0;
        while (!atrace[unit].empty() || !vtrace[unit].empty() || !ttrace[unit].empty() || !ftrace[unit].empty() || !ptrace[unit].empty())
        {
            pair<double, double> activity = atrace[unit].empty() ? missing: atrace[unit].front();
            pair<double, double> voltage = vtrace[unit].empty() ? missing: vtrace[unit].front();
            pair<double, double> temperature = ttrace[unit].empty() ? missing : ttrace[unit].front();
            pair<double, double> frequency = ftrace[unit].empty() ? missing : ftrace[unit].front();
            pair<double, double> power = ptrace[unit].empty() ? missing : ptrace[unit].front();
            
            // Expecting time in seconds
            double time = min({activity.first, voltage.first, temperature.first, frequency.first, power.first});
            DataPoint point(time - prev_time);
            if (activity.first == time)
            {
                point.activity = activity.second;
                atrace[unit].pop_front();
            }
            if (voltage.first == time)
            {
                point.vdd = voltage.second; // Expecting voltage in V
                vtrace[unit].pop_front();
            }
            if (temperature.first == time)
            {
                point.temperature = temperature.second; // Expecting temperature in K
                ttrace[unit].pop_front();
            }
            if (frequency.first == time)
            {
                point.frequency = frequency.second*1e6; // Expecting frequency in MHz, have to convert to Hz
                ftrace[unit].pop_front();
            }
            if (power.first == time)
            {
                point.power = power.second; // Expecting power in W
                ptrace[unit].pop_front();
            }
            trace[unit].push_back(point);

            prev_time = time;
        }
    }
    return trace;
}