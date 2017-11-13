#include "util.hh"

#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace oldspot
{

vector<string> split(const string& str, char delimiter)
{
    string token;
    istringstream stream(str);
    vector<string> tokens;
    while (getline(stream, token, delimiter))
        tokens.push_back(token);
    return tokens;
}

}