#ifndef UTILS_HPP
#define UTILS_HPP

#include <ostream>

struct Interval {
  double min;
  double max;
};

struct Point
{
  double x = 0;
  double y = 0;
  double z = 0;

  friend bool operator==(const Point& lhs, const Point& rhs) {
    /*if(std::abs(lhs.x-rhs.x) < std::numeric_limits<double>::epsilon() &&
     std::abs(lhs.y-rhs.y) < std::numeric_limits<double>::epsilon() &&
     std::abs(lhs.z-rhs.z) < std::numeric_limits<double>::epsilon()) {
    */
    if(std::abs(lhs.x-rhs.x) <= lhs.x/1E14 && std::abs(lhs.y-rhs.y) <= lhs.y/1E14 && std::abs(lhs.z-rhs.z) <= lhs.z/1E14)
      return true;

    return false; 
  }

  friend std::ostream& operator<<(std::ostream& os, const Point& obj) {
    os << '(' << obj.x << ", " << obj.y << ", " << obj.z << ')';
    return os;
  }
};
  
inline double absolute(double x) {
  if (x < 0)
    return -x;
  return x;
}

#endif //UTILS_HPP