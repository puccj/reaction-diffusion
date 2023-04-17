//#include "matrix.h"
#include "simulation.h"
#include <iostream>

int main() {
  /*
  double prova[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  Matrix m(4, 5, prova);
  m.print();
  
  double a = m[3][2];
  std::cout << a << '\n';

  a = -3;

  m[2][3] = a;

  m.print();
  */

  Simulation s(0,10,0,10, 0.1);
  s.setDu(1);
  s.setk2(11);
  s.setDv(0.04);
  s.setk1(7);
  
  s.lastU().print();
}