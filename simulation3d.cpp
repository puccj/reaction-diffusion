#include "simulation3d.h"
#include <random>
#include <fstream>

void Simulation3D::createFirstTensors(double* extremes /*, Func firstU, Func firstV*/) {
  int I = (extremes[1] - extremes[0]) / _h; 
  int J = (extremes[3] - extremes[2]) / _h; 
  int J = (extremes[5] - extremes[4]) / _h; 

  Tensor u(I,J,K);
  Tensor v(I,J,K);

  std::random_device gen;
  std::uniform_real_distribution<double> dist(-1,1);

  for (int i = 0; i < I; ++i)
    for (int j = 0; j < J; ++j) {
      for (int k = 0; j < K; ++k) {
        //default functions
        u[i][j][k] = 1+ 0.04*_k2*_k2 + 0.1*dist(gen);
        v[i][j][k] = 0.2*_k2 + 0.1*dist(gen);
      }
    }

  _u.push_back(std::move(u));
  _v.push_back(std::move(v));
}

Simulation3D::Simulation3D(double iMin, double iMax, double jMin, double jMax, double kMin, double kMax, double h, double k2)
    : _h{h}, _k2{k2} {
  double arr[6] = {iMin, iMax, jMin, jMax, kMin, kMax};
  createFirstTensores(arr);
}

Simulation3D::Simulation3D(Interval i, Interval j, Interval k, double h) : _h{h} {
  double arr[4] = {i.min, i.max, j.min, j.max, k.min, k.max};
  createFirstTensores(arr);
}

void Simulation3D::evolve(double dt, bool overwrite) {
  int size = _u.size(); // = _v.size()
  
  #define u _u[size-1]
  #define v _v[size-1]

  int I = u.i();  // = v.i()
  int J = u.j();  // = v.j()
  int K = u.k();  // = v.k()

  Tensor nextU(I,J,K);
  Tensor nextV(I,J,K);

  for (int i = 0; i < I; ++i) {
    for (int j = 0; j < J; ++j) {
      for (int k = 0; k < K; ++k) {
        double u_ijk = u[i][j][k];
        double v_ijk = v[i][j][k];
        double den = 1 + v_ijk*v_ijk;
        nextU[i][j][k] = u_ijk + (_Du*u.der2(i,j,k)/(_h*_h) + _k1 *(v_ijk - u_ijk*v_ijk /den)) *dt;
        nextV[i][j][k] = v_ijk + (_Dv*v.der2(i,j,k)/(_h*_h) + _k2 - v_ijk - 4*u_ijk*v_ijk /den) *dt;
      }
    }
  }

  if (overwrite) {
    _u[size-1] = std::move(nextU);
    _v[size-1] = std::move(nextV);
  }
  else {
    _u.push_back(std::move(nextU));
    _v.push_back(std::move(nextV));
  }
}

void Simulation3D::saveV(std::string filename){
  std::cout << "Saving data as '" << filename << "'... ";
  std::fstream fout(filename, std::ios::out);
  fout << _v.size() << '\n';
  fout << _v[0].rows() << ' ' << _v[0].cols() << '\n';

  for (auto m : _v) {
    fout << m << '\n';
  }
  fout.close();
  
  std::cout << "     Done!\n";
}
