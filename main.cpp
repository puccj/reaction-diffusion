#include "simulation3.h"
#include <iostream>

int main() {

  Simulation3D s({0,10},{0,10},{0,10},1);
  s.saveSurface("surface1.dat");
  s.saveMap("map1.dat");

  /*
  Simulation3D s("map05.dat", 1, 0.03, 5, 11);
  
  s.saveV("0.dat");
  std::cout << "Done 0/4\n";


  for (int i = 0; i < 20000; ++i) {
    s.evolve();
    // if (i % 50 == 0)
    //   s.saveV("evolved" + std::to_string(i) + ".dat");
  }
  s.saveV("2000.dat");
  std::cout << "Done 1/4\n";

  for (int i = 0; i < 120000; ++i)
    s.evolve();
  s.saveV("14000.dat");
  std::cout << "Done 2/4\n";

  for (int i = 0; i < 860000; ++i)
    s.evolve();
  s.saveV("100000.dat");
  std::cout << "Done 3/4\n";
  
  for (int i = 0; i < 1000000; ++i)
    s.evolve();
  s.saveV("200000.dat");
  std::cout << "Done 4/4\n";

  */
  std::cout << "Fine\n";
  return 0;
}