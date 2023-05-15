#include "simulation3.h"
#include <iostream>

int main() {
  Simulation3D s({0,10},{0,10},{0,10},0.5);
  s.saveMap();
  
  s.saveSurface();
  s.saveV();

  for (int i = 0; i < 1000; ++i) {
    s.evolve();
  }

  s.saveV("value2.dat");

  std::cout << "Fine\n";
  return 0;
}