#include "trace.hh"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "util.hh"

namespace oldspot
{

using namespace std;

/**
 * Push a string representation of a DataPoint onto a stream.
 */
ostream& operator<<(ostream& stream, const DataPoint& point)
{
    stream << point.time << ":{";
    stream << accumulate(next(point.data.begin()), point.data.end(),
                         point.data.begin()->first + ':' + to_string(point.data.begin()->second),
                         [](string a, const pair<string, double>& b){ return a + ',' + b.first + ':' + to_string(b.second); });
    return stream << "}";
}

/**
 * Parse a trace file to get a list of DataPoints.  Each trace file should be a table
 * delimited with the given delimiter (default is comma) where each row is a point in
 * the trace and each column is a quantity in the trace, with the first row containing
 * headers.  The first column must be the time at which the data point occurs.
 */
vector<DataPoint> parseTrace(const string fname, char delimiter)
{
    ifstream file(fname);
    if (!file)
    {
        cerr << fname << ": unable to open file" << endl;
        exit(1);
    }

    string line;
    vector<DataPoint> trace;
    
    // Parse unit names
    getline(file, line); // Get line containing unit names
    vector<string> quantities = split(line, delimiter);
    quantities.erase(quantities.begin());

    // Parse times (first column) and values
    double prev = 0;
    while (getline(file, line))
    {
        unordered_map<string, double> data;
        vector<string> values = split(line, delimiter);
        double time = stod(values[0]); // First column should be time
        for (size_t i = 0; i < quantities.size(); i++)
            data[quantities[i]] = stod(values[i + 1]);
        trace.push_back({time, time - prev, data});
        prev = time;
    }

    return trace;
}

} // namespace oldspot