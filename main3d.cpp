#include "simulation3.h"
#include <iostream>

int main() {
  /*
  Simulation3D s({0,11},{0,11},{0,11},0.1);
  s.saveSurface("surface01.dat");
  s.saveMap("map01.dat");
  */

  double Du = 1;  //fixed
  double k2 = 11; //fixed
  for (double Dv = 0.01; Dv <= 0.06; Dv += 0.01) {
    for (int k1 = 1; k1 <= 12; k1 += 2) {

      // double Dv = 0.06;
      // int k1 = 3;

      std::string name = std::to_string((int)Du) + ",0.0" + std::to_string((int)(Dv*100)) + ',' + std::to_string((int)k1) + ',' + std::to_string((int)k2) + " - ";

      Simulation3D s("map01.dat", Du, Dv, k1, k2);
      
      // s.saveV("v - 0.dat");
      // s.saveU("u - 0.dat");

      s.saveU("data/3D/" + name + "v - 0.dat");
      s.saveV("data/3D/" + name + "u - 0.dat");

      for (int step = 100; step < 2002; step+=100) {
        for (int i = 0; i < 100; ++i) {
          s.evolve();
        }
        s.saveU("data/3D/" + name + "v - " + std::to_string(step) + ".dat");
        s.saveV("data/3D/" + name + "u - " + std::to_string(step) + ".dat");
      }      

      std::cout << "End\n";
    }
  }

  std::cout << "End All\n";
  return 0;
}