#ifndef SIMULATION_H
#define SIMULATION_H

#include "matrix.h"
#include <iostream>

struct Interval {
  double min;
  double max;
};

class Simulation {
 private:
  double _Du;
  double _Dv;
  double _k1;
  double _k2;
  double _h;  //Uniform mesh size
  //contains the matrix of u and v at each time step
  std::vector<Matrix> _u;
  std::vector<Matrix> _v;

  void createFirstMatrixes(double* extremes);
  /*
  template<typename Func>
  void funcToMatrix(Func initialU, Func initialV, double* extremes) {
    int rows = (extremes[3] - extremes[1]) / _h; //yMax - yMin
    int cols = (extremes[2] - extremes[0]) / _h; //xMax - xMin
    Matrix u(rows, cols);
    Matrix v(rows, cols);

    for (int i = 0; i < rows; ++i)
      for (int j = 0; j < cols; ++j) {
        u[i][j] = initialU(i,j);
        v[i][j] = initialV(i,j);
      }

    _u.push_back(u);
    _v.push_back(v);
  }
  */

 public:
  Simulation(Interval x, Interval y, /*Func firstU, Func firstV,*/ double h = 1);
  Simulation(double xMin = 0, double xMax = 10, double yMin = 0, double yMax = 10, double h = 1);
  Simulation(Matrix& initialU, Matrix& initialV, double h = 1) : _u{initialU}, _v{initialV}, _h{h} {};

  /*
  template<typename Func>
  Simulation(Func& initialU, Func& initialV, double xMin, double xMax, double yMin, double yMax, double h = 1) {
    double arr[4] = {xMin, xMax, yMin, yMax};
    funcToMatrix(initialU, initialV, arr);
  }

  template<typename Func>
  Simulation(Func& initialU, Func& initialV, Interval x, Interval y, double h = 1) {
    double arr[4] = {x.min, x.max, y.min, y.max};
    funcToMatrix(initialU, initialV, arr);
  }
  */

  void setDu(double val) {_Du = val; };
  void setDv(double val) {_Dv = val; };
  void setk1(double val) {_k1 = val; };
  void setk2(double val) {_k2 = val; };

  Matrix lastU() {return _u[_u.size()-1]; };
  Matrix lastV() {return _v[_v.size()-1]; };

  //calculate the next step of the Matrix, given the time interval dT
  void evolve(double dt);
};

#endif //SIMULATION_H