#ifndef SIMULATION_H
#define SIMULATION_H

#include "matrix.h"
#include "utils.hpp"
#include <iostream>
#include <vector>

class Simulation {
 private:
  Matrix _u; //concentration of the inhibitor
  Matrix _v; //concentration of the activator
  double _h; //Uniform mesh size
  double _Du;
  double _Dv;
  double _k1;
  double _k2;
  
 public:
  Simulation(double xMin = 0, double xMax = 10, double yMin = 0, double yMax = 10, double h = 1, double k2 = 11);
  Simulation(Interval x, Interval y, /*Func firstU, Func firstV,*/ double h = 1, double k2 = 11);
  Simulation(Matrix& initialU, Matrix& initialV, double h = 1) : _u{initialU}, _v{initialV}, _h{h} {};

  void setDu(double val) {_Du = val; };
  void setDv(double val) {_Dv = val; };
  void setk1(double val) {_k1 = val; };
  void setk2(double val);

  /*
  Matrix lastU() {return _u[_u.size()-1]; };
  Matrix lastV() {return _v[_v.size()-1]; };
  */

  //calculate the next step of the Matrix, given the time interval dT
  void evolve(double dt);

  //save activator concentration in a file
  void saveV(std::string filename = "simulation.dat");

 private:
  void fillUV();
};

#endif //SIMULATION_H