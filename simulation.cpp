#include "simulation.h"
#include <random>

void Simulation::createFirstMatrixes(double* extremes /*, Func firstU, Func firstV*/) {
  int rows = (extremes[3] - extremes[2]) / _h; //yMax - yMin
  int cols = (extremes[1] - extremes[0]) / _h; //xMax - xMin

  Matrix u(rows, cols);
  Matrix v(rows, cols);

  std::random_device gen;
  std::uniform_real_distribution<double> dist(-1,1);

  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j) {
      //default functions
      u[i][j] = 1+ 0.04*_k2*_k2 + 0.1*dist(gen);
      v[i][j] = 0.2*_k2 + 0.1*dist(gen);
    }

  _u.push_back(u);
  _v.push_back(v);
}

Simulation::Simulation(double xMin, double xMax, double yMin, double yMax, double h) : _h{h} {
  double arr[4] = { xMin, xMax, yMin, yMax};
  createFirstMatrixes(arr);
}

Simulation::Simulation(Interval x, Interval y, double h) : _h{h} {
  double arr[4] = {x.min, x.max, y.min, y.max};
  createFirstMatrixes(arr);
}

void Simulation::evolve(double dt) {
  int size = _u.size();   // = _v.size()
  Matrix u = _u[size-1];
  Matrix v = _v[size-1];

  int rows = u.rows();  // = lastV.rows()
  int cols = u.cols();  // = lastV.cols()
  Matrix nextU(u);
  Matrix nextV(v);

  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j) {
      double u_ij = u[i][j];
      double v_ij = v[i][j];
      double den = 1 + v_ij*v_ij;
      nextU[i][j] = u_ij + (_Du*u.der2(i,j)/(_h*_h) + _k1*(v_ij - u_ij*v_ij /den)) *dt;
      nextV[i][j] = v_ij + (_Dv*v.der2(i,j)/(_h*_h) + _k2 - v_ij - 4*u_ij*v_ij /den) *dt;
    }
}