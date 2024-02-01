#include "simulation3.h"
#include <iostream>

int main() {

  /*
  Simulation3D s({0,11},{0,11},{0,11},0.1);
  s.saveSurface("surface01.dat");
  s.saveMap("map01.dat");
  */ 

  Simulation3D s("map01.dat", 1, 0.03, 5, 11);
  
  s.saveV("v - 0.dat");
  s.saveU("u - 0.dat");

  for (int step = 1; step < 2002; ++step) {
    s.evolve();
    s.saveU("1,0.03,5,11 - v - " + std::to_string(step) + ".dat");
    s.saveV("1,0.03,5,11 - u - " + std::to_string(step) + ".dat");
  }

  /*
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

  std::cout << "End\n";
  return 0;
}