#ifndef MESH_H
#define MESH_H
#include <ostream>

//To do: template<class T>
class Mesh {
  int _i;
  int _j;
  int _k;
  int _h; //uniform mesh size
  double* _data;

  struct Vec
  {
    int k;
    double* data = nullptr;
    double operator[](int index) const { return data[index]; };
    double& operator[](int index) { return data[index]; };
  };

  struct Mat
  {
    int j;
    int k;
    double* data = nullptr;
    Vec operator[](int index) const { return Vec{k, data + k*index}; };
    //Vec& operator[](int index) { return Vec{k, data + k*index}; };
  };

 public:
  Mesh(int i, int j, int k, double data = double{});
  Mesh(int i, int j, int k, double *data);
  Mesh(const Mesh &src);            //copy constructor
  Mesh(Mesh &&src);			            //move constructor
  Mesh& operator=(const Mesh &src); //copy assignment
  Mesh& operator=(Mesh &&src);	    //move assignment
  ~Mesh();                          //destructor

  int i() { return _i; };
  int j() { return _j; };
  int k() { return _k; };
  int h() { return _h; };
  Mat operator[](int index) const { return Mat{_j, _k, _data + _j*_k*index}; };

  //print the sensor through slices
  void print();
  
  //return the discretized second derivative at point i,j,k
  double der2(int i, int j, int k);
  //return a Mesh whose values are the discretized second derivative
  Mesh der2();
  
  //overload of (some) operators
  /*
  Mesh& operator+=(const Mesh& rhs);
  Mesh& operator-=(const Mesh& rhs);
  Mesh& operator*=(const double& rhs);
  friend Mesh operator*(Mesh lhs, double const& rhs);
  friend Mesh operator*(double const& lhs, Mesh rhs);
  */
  friend std::ostream& operator<<(std::ostream& os, const Mesh& obj);
};

#endif //MESH_H