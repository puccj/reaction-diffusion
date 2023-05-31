#include "surface.h"
#include <limits>
#include <iostream>
#include <iomanip>

Surface::Surface(int nPoints, Point data) :
                _nPoints{nPoints},
                _data{new Point[nPoints]} {
  std::fill(_data, _data + nPoints, data);
}

Surface::Surface(int nPoints, Point *data) :
                _nPoints{nPoints},
                _data{new Point[nPoints]} {
  std::copy(data, data + nPoints, _data);
}

Surface::Surface(const Surface &src) :
                _nPoints{src._nPoints},
                _data{new Point[src._nPoints]} {
  
  std::copy(src._data, src._data + src._nPoints, _data);

}

Surface::Surface(Surface &&src) : 
              _nPoints{src._nPoints},
              _data{src._data} {
  src._nPoints = 0;
  src._data = nullptr;
}

Surface &Surface::operator=(const Surface &src) {
  // Guard self assignment
  if (this == &src)
    return *this;

  if (_nPoints != src._nPoints) {           // resource in *this cannot be reused
    Point* temp = new Point[src._nPoints];  // allocate resource, if throws, do nothing
    delete[] _data;                         // release resource in *this
    _data = temp;
    _nPoints = src._nPoints;
  } 

  std::copy(src._data, src._data + src._nPoints, _data);
  return *this;
}

Surface &Surface::operator=(Surface &&src) {
  // Guard self assignment
  if (this == &src)
    return *this;

  delete[] _data;       // release resource in *this
  _data = src._data;
  _nPoints = src._nPoints;
  src._data = nullptr;  // leave src in valid state
  src._nPoints = 0;

  return *this;
}

Surface::~Surface() {
  delete[] _data;
}

std::ostream &operator<<(std::ostream &os, const Surface &obj)
{
  //format: 3 cols of value: x y z
  for (int i = 0; i < obj._nPoints; ++i) {
    os << obj._data[i].x << ' ' << obj._data[i].y << ' ' << obj._data[i].z << '\n';
  }
  return os;
}