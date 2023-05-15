#ifndef SIMULATION_H
#define SIMULATION_H

#include "matrix.h"
#include <iostream>
#include <vector>

struct Interval {
  double min;
  double max;
};

class Simulation {
 private:
  std::vector<Matrix> _u; //u(x,t): concentration of the inhibitor
  std::vector<Matrix> _v; //v(x,t): concentration of the activator
  double _h;              //Uniform mesh size
  double _Du;
  double _Dv;
  double _k1;
  double _k2;

  void createFirstMatrixes(double* extremes);
  
 public:
  Simulation(double xMin = 0, double xMax = 10, double yMin = 0, double yMax = 10, double h = 1, double k2 = 11);
  Simulation(Interval x, Interval y, /*Func firstU, Func firstV,*/ double h = 1);
  Simulation(Matrix& initialU, Matrix& initialV, double h = 1) : _u{initialU}, _v{initialV}, _h{h} {};

  void setDu(double val) {_Du = val; };
  void setDv(double val) {_Dv = val; };
  void setk1(double val) {_k1 = val; };
  void setk2(double val) {_k2 = val; };

  Matrix lastU() {return _u[_u.size()-1]; };
  Matrix lastV() {return _v[_v.size()-1]; };

  //calculate the next step of the Matrix, given the time interval dT
  void evolve(double dt, bool overwrite = false);

  //save activator concentration in a file
  void saveV(std::string filename = "simulation.dat");
};

#endif //SIMULATION_H