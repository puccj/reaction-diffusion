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
  return _data[_cols*(i+1) +j] + _data[_cols*(i-1) +j] + 
    _data[_cols*i +j+1] + _data[_cols*i +j-1] - 6*_data[_cols*i +j];
}

Matrix Matrix::der2() {
  Matrix derivative(_rows, _cols);

  //To do: manage boundaries conditions
  for (int i = 1; i < _rows-1; ++i) {
    for (int j = 1; j < _cols-1; ++j) {
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