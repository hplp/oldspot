#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace oldspot
{

template<typename T>
inline T linterp(const T& x, const std::pair<T, T>& s, const std::pair<T, T>& f)
{
    return s.second + (f.second - s.second)*(x - s.first)/(f.first - s.first);
}

std::vector<std::string> split(const std::string& str, char delimiter);

template<typename T>
void print_table(const std::vector<std::string>& rows, const std::vector<std::string>& cols, const std::unordered_map<std::string, std::unordered_map<std::string, T>>& data)
{
    using namespace std;

    unordered_map<string, size_t> col_widths;
    size_t first_col_width = 0;
    for (const string& header: cols)
        col_widths[header] = header.length();
    for (const string& header: rows)
        first_col_width = max(first_col_width, header.length());
    for (size_t i = 0; i < rows.size(); i++)
        for (size_t j = 0; j < cols.size(); j++)
            col_widths[cols[j]] = max(col_widths[cols[j]], to_string(data.at(rows[i]).at(cols[j])).length());

    // Print column headers
    for (size_t k = 0; k < first_col_width; k++)
        cout << ' ';
    cout << " | ";
    for (size_t j = 0; j < cols.size(); j++)
        cout << left << setw(col_widths[cols[j]]) << cols[j] << " | ";
    cout << endl;
    // Print row data
    for (size_t i = 0; i < rows.size(); i++)
    {
        cout << left << setw(first_col_width) << rows[i] << " | ";
        for (size_t j = 0; j < cols.size(); j++)
            cout << right << setw(col_widths[cols[j]]) << to_string(data.at(rows[i]).at(cols[j])) << " | ";
        cout << endl;
    }
}

} // namespace oldspot