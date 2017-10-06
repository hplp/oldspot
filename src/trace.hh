#pragma once

#include <cstdlib>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

struct DataPoint
{
    double time;
    std::map<std::string, double> data;

    friend std::ostream& operator<<(std::ostream& stream, const DataPoint& point);
};

std::vector<DataPoint> parseTrace(const std::string fname, char delimiter=',');