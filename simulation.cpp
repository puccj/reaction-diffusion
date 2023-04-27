#include "simulation.h"
#include <random>
#include <fstream>

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

Simulation::Simulation(double xMin, double xMax, double yMin, double yMax, double h, double k2) : _h{h}, _k2{k2} {
  double arr[4] = { xMin, xMax, yMin, yMax};
  createFirstMatrixes(arr);
}

Simulation::Simulation(Interval x, Interval y, double h) : _h{h} {
  double arr[4] = {x.min, x.max, y.min, y.max};
  createFirstMatrixes(arr);
}

void Simulation::evolve(double dt, bool keep) {
  int size = _u.size(); // = _v.size()
  
  std::cout << "Debug: matrix u" << _u[size-1] << '\n';
  std::cout << "Debug: matrix v" << _v[size-1] << '\n';
  
  Matrix u = _u[size-1];
  Matrix v = _v[size-1];

  int rows = u.rows();  // = v.rows()
  int cols = u.cols();  // = v.cols()

  Matrix nextU(rows, cols);
  Matrix nextV(rows, cols);

  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j) {
      double u_ij = u[i][j];
      double v_ij = v[i][j];
      double den = 1 + v_ij*v_ij;
      nextU[i][j] = u_ij + (_Du*u.der2(i,j)/(_h*_h) + _k1*(v_ij - u_ij*v_ij /den)) *dt;
      nextV[i][j] = v_ij + (_Dv*v.der2(i,j)/(_h*_h) + _k2 - v_ij - 4*u_ij*v_ij /den) *dt;
    }

  if (keep) {
    _u.push_back(nextU);
    _v.push_back(nextV);
  }
  else {
    _u[size-1] = nextU;
    _v[size-1] = nextV;
  }

  std::cout << "Debug: matrix u" << _u[size-1] << '\n';
  std::cout << "Debug: matrix v" << _v[size-1] << '\n';
}

void Simulation::saveV(std::string filename){
  std::cout << "Saving data...";
  std::fstream fout(filename, std::ios::out);
  fout << _v.size() << '\n';
  fout << _v[0].rows() << ' ' << _v[0].cols() << '\n';

  for (auto m : _v) {
    fout << m << '\n';
  }
  fout.close();
  
  std::cout << "     Done!\n";
}
