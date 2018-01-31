#pragma once

#include <cstdlib>
#include <limits>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace oldspot
{

/**
 * Data point in an activity trace for a unit.  It contains a time at which it
 * occurs, the duration of the segment, and a map of quantities onto their values
 * for that segment (i.e. temperature, voltage, frequency, etc., depending on
 * which quantities are needed to compute reliability).
 */
struct DataPoint
{
    double time;
    double duration;
    std::unordered_map<std::string, double> data;

    friend std::ostream& operator<<(std::ostream& stream, const DataPoint& point);
};

std::vector<DataPoint> parseTrace(const std::string fname, char delimiter=',');

} // namespace oldspot