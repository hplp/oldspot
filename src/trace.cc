#include "trace.hh"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace oldspot
{

using namespace std;

ostream& operator<<(ostream& stream, const DataPoint& point)
{
    stream << point.time << ":{";
    stream << accumulate(next(point.data.begin()), point.data.end(),
                         point.data.begin()->first + ':' + to_string(point.data.begin()->second),
                         [](string a, const pair<string, double>& b){ return a + ',' + b.first + ':' + to_string(b.second); });
    return stream << "}";
}

vector<DataPoint> parseTrace(const string fname, char delimiter)
{
    ifstream file(fname);
    if (!file)
    {
        cerr << fname << ": unable to open file" << endl;
        exit(1);
    }

    string line, token;
    istringstream stream;
    vector<string> quantities;
    vector<DataPoint> trace;
    
    // Parse unit names
    getline(file, line); // Get line containing unit names
    stream.str(line);
    stream.clear();
    getline(stream, token, delimiter); // Ignore header of time column
    while (getline(stream, token, delimiter))
        quantities.push_back(token);

    // Parse times (first column) and values
    double prev = 0;
    while (getline(file, line))
    {
        stream.str(line);
        stream.clear();
        getline(stream, token, delimiter);
        double time = stod(token); // First column should be time
        map<string, double> data;
        for (const string& quantity: quantities)
        {
            getline(stream, token, delimiter);
            data[quantity] = stod(token);
        }
        trace.push_back({time, time - prev, data});
        prev = time;
    }

    return trace;
}

} // namespace oldspot