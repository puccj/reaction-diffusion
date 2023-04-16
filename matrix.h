#ifndef MATRIX_H
#define MATRIX_H
#include <vector>

//TO DO: bound checking
  
//To do: template<class T>
class Matrix {
  int _rows;
  int _cols;
  double* _data;

  struct Vec
  {
    int lenght;
    double* data = nullptr;
    double operator[](int index) const { return data[index]; };
    double& operator[](int index) { return data[index]; };
  };

 public:
  Matrix(int rows, int cols, double data = double{});
  Matrix(int rows, int cols, double *data);
  ~Matrix();

  int rows() { return _rows; };
  int cols() { return _cols; };
  Vec operator[](int index) const { return Vec{_cols, _data + _cols*index}; };
  //Vec& operator[](int index) { return Vec{_cols, _data + _cols*index}; };

  //print the matrix
  void print();
  
  //return the discretized second derivative at point i,j
  double der2(int i, int j);
  //return a matrix whose values are the discretized second derivative
  Matrix der2();
  
  //Overload of (some) operators
  Matrix& operator+=(const Matrix& rhs);
  Matrix& operator-=(const Matrix& rhs);
  Matrix& operator*=(const double& rhs);
  friend Matrix operator*(Matrix lhs, double const& rhs);
  friend Matrix operator*(double const& lhs, Matrix rhs);
};

#endif //MATRIX_H