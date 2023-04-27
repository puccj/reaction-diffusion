#include "matrix.h"
#include <iostream>
#include <algorithm>

Matrix::Matrix(int rows, int cols, double data) :
              _rows{rows}, 
              _cols{cols}, 
              _data{new double[rows*cols]} {

  std::fill(_data, _data + rows*cols, data);
}

Matrix::Matrix(int rows, int cols, double *data) :
              _rows{rows}, 
              _cols{cols}, 
              _data{new double[rows*cols]} {
  
  std::copy(data, data + rows*cols, _data);
}

Matrix::~Matrix() {    
  delete[] _data;
}

void Matrix::print() {
  for (int i = 0; i < _rows; ++i) {
    for (int j = 0; j < _cols; ++j) {
      std::cout << _data[_cols*i+j] << "   ";
    }
    std::cout << '\n';
  }
}

double Matrix::der2(int i, int j) {
  if (i < 0 || i >= _rows || j < 0 || j > _cols) {
    std::cout << "Debug: i = " << i << " - j = " << j << " - rows = " << _rows << " - cols = " << _cols << '\n';
    std::cerr << "Error in der2 function: indexes not valid\n";
    return 0;
  }

  //Periodic boundary conditions
  double down, up, right, left;
  if (i == _rows-1)
    down = _data[j];
  else
    down = _data[_cols*(i+1) +j];

  if (i == 0)
    up = _data[_cols*(_rows-1) + j];
  else
    up = _data[_cols*(i-1) +j];

  if (j == _cols-1)
    right = _data[_cols*i];
  else
    right = _data[_cols*i +j+1];

  if (j == 0)
    left = _data[_cols*i + _cols-1];
  else
    left = _data[_cols*i + j-1];
    
  return down + up + right + left - 4*_data[_cols*i +j];
}

Matrix Matrix::der2() {
  Matrix derivative(_rows, _cols);

  for (int i = 0; i < _rows; ++i) {
    for (int j = 0; j < _cols; ++j) {
      derivative[i][j] = der2(i,j);
    }
  }

  return derivative;
}

Matrix& Matrix::operator+=(const Matrix& rhs) {
  for (int i = 0; i < _rows*_cols; ++i)
    _data[i] += rhs._data[i];
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs) {
  for (int i = 0; i < _rows*_cols; ++i)
    _data[i] -= rhs._data[i];
  return *this;
}

Matrix& Matrix::operator*=(const double& rhs) {
  for (int i = 0; i < _rows*_cols; ++i)
    _data[i] *= rhs;
  return *this;
}

Matrix operator*(Matrix lhs, const double& rhs) {
  lhs *= rhs;
  return lhs;
}

Matrix operator*(const double& lhs, Matrix rhs) {
  rhs *= lhs;
  return rhs;
}

std::ostream& operator<<(std::ostream& os, const Matrix& obj)
{
  for (int i = 0; i < obj._rows; ++i) {
    for (int j = 0; j < obj._cols; ++j) {
      os << obj._data[obj._cols*i+j] << " ";
    }
    os << '\n';
  }

  return os;
}