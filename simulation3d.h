#ifndef SIMULATION3D_H
#define SIMULATION3D_H

#include "tensor.h"
#include <iostream>
#include <vector>

struct Interval {
  double min;
  double max;
};

class Simulation3D {
 private:
  std::vector<Tensor> _u; //u(x,t): concentration of the inhibitor
  std::vector<Tensor> _v; //v(x,t): concentratiom of the activator
  double _h;              //Uniform mesh size
  double _Du;
  double _Dv;
  double _k1;
  double _k2;

  void createFirstTensors(double* extremes);
  
 public:
  Simulation3D(double iMin = 0, double iMax = 10, double jMin = 0, double jMax = 10, double kMin, double kMax, double h = 1, double k2 = 11);
  Simulation3D(Interval i, Interval j, Interval k, /*Func firstU, Func firstV,*/ double h = 1);
  Simulation3D(Tensor& initialU, Tensor& initialV, double h = 1) : _u{initialU}, _v{initialV}, _h{h} {};

  void setDu(double val) {_Du = val; };
  void setDv(double val) {_Dv = val; };
  void setk1(double val) {_k1 = val; };
  void setk2(double val) {_k2 = val; };

  Tensor lastU() {return _u[_u.size()-1]; };
  Tensor lastV() {return _v[_v.size()-1]; };

  //calculate the next step of the Tensor, given the time interval dT
  void evolve(double dt, bool overwrite = false);

  //save activator concentration in a file
  void saveV(std::string filename = "simulation3D.dat");
};

#endif //SIMULATION3D_H