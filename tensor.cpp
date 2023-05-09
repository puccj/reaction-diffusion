#include "tensor.h"
#include <iostream>
#include <algorithm>

Tensor::Tensor(int i, int j, int k, double data) :
              _i{i}, 
              _j{j},
              _k{k},
              _data{new double[i*j*k]} {

  std::fill(_data, _data + i*j*k, data);
}

Tensor::Tensor(int i, int j, int k, double *data) :
              _i{i}, 
              _j{j},
              _k{k},
              _data{new double[i*j*k]} {
  
  std::copy(data, data + i*j*k, _data);
}

Tensor::Tensor(const Tensor &src) : 
              _i{src._i},
              _j{src._j},
              _k{src._k},
              _data{new double[_i*_j*_k]} {
  
  std::copy(src._data, src._data + _i*_j*_k, _data);
}

Tensor::Tensor(Tensor &&src) : 
              _i{src._i}, 
              _j{src._j}, 
              _k{src._k},
              _data{src._data} {
  src._i = 0;
  src._j = 0;
  src._k = 0;
  src._data = nullptr;
}

Tensor &Tensor::operator=(const Tensor &other) {
  // Guard self assignment
  if (this == &other)
    return *this;

  size_t thisSize = _i * _j * _k;
  size_t otherSize = other._i * other._j * other._k;

  if (thisSize != otherSize) {            // resource in *this cannot be reused
    double* temp = new double[otherSize]; // allocate resource, if throws, do nothing
    delete[] _data;                       // release resource in *this
    _data = temp;
    _i = other._i;
    _j = other._j;
    _k = other._k;
  } 

  std::copy(other._data, other._data + otherSize, _data);
  return *this;
}

Tensor &Tensor::operator=(Tensor &&other) {
  // Guard self assignment
  if (this == &other)
    return *this;

  delete[] _data;         // release resource in *this
  _data = other._data;
  _i = other._i;
  _j = other._j;
  _k = other._k;
  other._data = nullptr;  // leave other in valid state
  other._i = 0;
  other._j = 0;
  other._k = 0;

  return *this;
}

Tensor::~Tensor() {
  delete[] _data;
}

void Tensor::print() {
  for (int i = 0; i < _i; ++i) {
    std::cout << "i = " << i << ":\n";

    for (int j = 0; j < _j; ++j) {
      for (int k = 0; k < _k; ++k) {
        std::cout << _data[i*_j*_k + j*_k +k] << "   ";
      }
    }
    std::cout << '\n';
  }
}

double Tensor::der2(int i, int j, int k) {
  if (i < 0 || i >= _i || j < 0 || j >= _j || k < 0 || k >= _k) {
    std::cerr << "Error in der2 function: indexes not valid\n";
    return 0;
  }

  //Periodic boundary conditions
  double back, front, down, up, right, left;
  if (i == _i-1)
    back = _data[j*_k + k];
  else
    back = _data[(i+1)*_j*_k + j*_k + k];

  if (i == 0)
    front = _data[(_i-1)*_j*_k + j*_k + k];
  else
    front = _data[(i-1)*_j*_k + j*_k + k];

  if (j == _j-1)
    down = _data[i*_j*_k + k];
  else
    down = _data[i*_j*_k + (j+1)*_k + k];

  if (j == 0)
    up = _data[i*_j*_k + (_j-1)*_k + k];
  else
    up = _data[i*_j*_k + (j-1)*_k + k];

  if (k == _k-1)
    right = _data[i*_j*_k + j*_k];
  else
    right = _data[i*_j*_k + j*_k + k+1];

  if (k == 0)
    left = _data[i*_j*_k + j*_k + _k-1];
  else
    left = _data[i*_j*_k + j*_k + k-1];
    
  return down + up + right + left - 6*_data[i*_j*_k + j*_k + k];
}

Tensor Tensor::der2() {
  Tensor derivative(_i, _j, _k);

  for (int i = 0; i < _i; ++i) {
    for (int j = 0; j < _j; ++j) {
      for (int k = 0; k < _k; ++k) {
        derivative[i][j][k] = der2(i,j,k);
      }
    }
  }

  return derivative;
}
/*
Tensor& Tensor::operator+=(const Tensor& rhs) {
  for (int i = 0; i < _i*_j; ++i)
    _data[i] += rhs._data[i];
  return *this;
}

Tensor& Tensor::operator-=(const Tensor& rhs) {
  for (int i = 0; i < _i*_j; ++i)
    _data[i] -= rhs._data[i];
  return *this;
}

Tensor& Tensor::operator*=(const double& rhs) {
  for (int i = 0; i < _i*_j; ++i)
    _data[i] *= rhs;
  return *this;
}

Tensor operator*(Tensor lhs, const double& rhs) {
  lhs *= rhs;
  return lhs;
}

Tensor operator*(const double& lhs, Tensor rhs) {
  rhs *= lhs;
  return rhs;
}
*/

std::ostream& operator<<(std::ostream& os, const Tensor& obj)
{
  for (int i = 0; i < obj._i; ++i) {
    os << "i = " << i << ":\n";
    for (int j = 0; j < obj._j; ++j) {
      for (int k = 0; k < obj._k; ++k) {
        os << obj._data[i*obj._j*obj._k + j*obj._k + k] << " ";
      }
      os << '\n';
    }
    os << '\n';
  }

  return os;
}