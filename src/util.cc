#include "util.hh"

#include <sstream>
#include <string>
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

}