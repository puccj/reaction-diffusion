#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include <ostream>

struct Point
{
  double x = 0;
  double y = 0;
  double z = 0;
  friend bool operator==(const Point& lhs, const Point& rhs);
  friend std::ostream& operator<<(std::ostream& os, const Point& obj);
};

//to do: template class T
class Surface {
  int _nPoints;
  Point* _data;
  
 public:
  Surface(int nPoints = 0, Point data = Point{});
  Surface(int nPoints, Point* data);
  Surface(const Surface &src);            //copy constructor
  Surface(Surface &&src);			            //move constructor
  Surface& operator=(const Surface &src); //copy assignment
  Surface& operator=(Surface &&src);	    //move assignment
  ~Surface();                             //destructor

  int nPoints() { return _nPoints; };
  Point operator[](int index) const { return _data[index]; };

  friend std::ostream& operator<<(std::ostream& os, const Surface& obj);
};

#endif  //SURFACE_H