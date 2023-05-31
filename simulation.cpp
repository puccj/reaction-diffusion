#include "simulation.h"
#include <random>
#include <fstream>

Simulation::Simulation(double xMin, double xMax, double yMin, double yMax, double h, double k2) 
  : _h{h}, 
    _k2{k2}
{
  _u = std::move(Matrix((int)(yMax-yMin) /_h, (int)(xMax-xMin) /_h)),
  _v = std::move(Matrix((int)(yMax-yMin) /_h, (int)(xMax-xMin) /_h)),
  fillUV();
}

Simulation::Simulation(Interval x, Interval y, double h, double k2) 
  : _h{h},
    _k2{k2} 
{
  _u = std::move(Matrix((int)(y.max-y.min) /_h, (int) (x.max-x.min) /_h)),
  _v = std::move(Matrix((int)(y.max-y.min) /_h, (int) (x.max-x.min) /_h)),
  fillUV();
}

void Simulation::setk2(double val) {
  _k2 = val;
  fillUV();
}

void Simulation::evolve(double dt) {
  int rows = _u.rows();  // == v.rows()
  int cols = _u.cols();  // == v.cols()

  Matrix nextU(rows, cols);
  Matrix nextV(rows, cols);

  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j) {
      double u_ij = _u[i][j];
      double v_ij = _v[i][j];
      double den = 1 + v_ij*v_ij;
      nextU[i][j] = u_ij + (_Du*_u.der2(i,j)/(_h*_h) + _k1*(v_ij - u_ij*v_ij /den)) *dt;
      nextV[i][j] = v_ij + (_Dv*_v.der2(i,j)/(_h*_h) + _k2 - v_ij - 4*u_ij*v_ij /den) *dt;
    }

    _u = std::move(nextU);
    _v = std::move(nextV);
}

void Simulation::saveV(std::string filename) {
  std::fstream fout(filename, std::ios::out);
  //fout << _v.rows() << ' ' << _v.cols() << '\n';
  fout << _v << '\n';
  fout.close();
}


//private

void Simulation::fillUV() {
  int rows = _u.rows();   //== _v.rows
  int cols = _u.cols();

  std::cout << "Debug: rows = " << rows << '\n';
  std::cout << "Debug: cols = " << cols << '\n';

  std::random_device gen;
  std::uniform_real_distribution<double> dist(-1,1);
  
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      //default functions
      _u[i][j] = 1+ 0.04*_k2*_k2 + 0.1*dist(gen);
      _v[i][j] = 0.2*_k2 + 0.1*dist(gen);
    }
  }
}
