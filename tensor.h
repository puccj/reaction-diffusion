#ifndef TENSOR_H
#define TENSOR_H
#include <ostream>

//To do: template<class T>
class Tensor {
  int _i;
  int _j;
  int _k;
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
  Tensor(int i, int j, int k, double data = double{});
  Tensor(int i, int j, int k, double *data);
  Tensor(const Tensor &src);            //copy constructor
  Tensor(Tensor &&src);			            //move constructor
  Tensor& operator=(const Tensor &src); //copy assignment
  Tensor& operator=(Tensor &&src);	    //move assignment
  ~Tensor();                            //destructor

  int i() { return _i; };
  int j() { return _j; };
  int k() { return _k; };
  Mat operator[](int index) const { return Mat{_j, _k, _data + _j*_k*index}; };

  //print the sensor through slices
  void print();
  
  //return the discretized second derivative at point i,j,k
  double der2(int i, int j, int k);
  //return a Tensor whose values are the discretized second derivative
  Tensor der2();
  
  //overload of (some) operators
  /*
  Tensor& operator+=(const Tensor& rhs);
  Tensor& operator-=(const Tensor& rhs);
  Tensor& operator*=(const double& rhs);
  friend Tensor operator*(Tensor lhs, double const& rhs);
  friend Tensor operator*(double const& lhs, Tensor rhs);
  */
  friend std::ostream& operator<<(std::ostream& os, const Tensor& obj);
};

#endif //TENSOR_H