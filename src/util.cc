#include "util.hh"

#include <cstdio>
#include <cstdarg>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

using namespace std;

namespace oldspot
{

/**
 * Split a string into tokens with the given single-character delimiter.
 */
vector<string>
split(const string& str, char delimiter)
{
    if (str.empty())
        return {""};

    string token;
    istringstream stream(str);
    vector<string> tokens;
    while (getline(stream, token, delimiter))
        tokens.push_back(token);
    return tokens;
}

int warn(const char* format, ...)
{
    va_list args1;
    va_start(args1, format);
    va_list args2;
    va_copy(args2, args1);
    vector<char> buf(vsnprintf(NULL, 0, format, args1) + 1);
    va_end(args1);
    vsnprintf(buf.data(), buf.size(), format, args2);
    va_end(args2);

    static unordered_set<string> warned;
    string str(buf.begin(), prev(buf.end()));
    if (warned.count(str) == 0)
    {
        warned.insert(str);
        cerr << "warning: " << str;
    }

    return str.size();
}

}