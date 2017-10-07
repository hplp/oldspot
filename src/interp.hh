#include <utility>

template<typename T>
inline T linterp(const T& x, const std::pair<T, T>& s, const std::pair<T, T>& f)
{
    return s.second + (f.second - s.second)*(x - s.first)/(f.first - s.first);
}