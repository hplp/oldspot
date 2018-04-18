#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "unit.hh"

namespace oldspot
{

/**
 * Linearly interpolate between two points, s (start) and f (end).
 */
template<typename T> inline T
linterp(const T& x, const std::pair<T, T>& s, const std::pair<T, T>& f)
{
    return s.second + (f.second - s.second)*(x - s.first)/(f.first - s.first);
}

std::vector<std::string> split(const std::string& str, char delimiter);

int warn(const char* format, ...);

template<typename Rows> void
writecsv(const std::string& filename, const std::vector<std::shared_ptr<Unit>>& units, Rows&& rows)
{
    using namespace std;

    ofstream data(filename);
    if (data)
    {
        for (const auto& row: rows)
            data << ',' << row.first;
        data << endl;
        for (const shared_ptr<Unit>& unit: units)
        {
            data << unit->name;
            for (const auto& row: rows)
                data << ',' << row.second(unit);
            data << endl;
        }
    }
    else
        cerr << "error: could not write to " << filename << endl;
}

} // namespace oldspot