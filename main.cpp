#include "simulation3.h"
#include <iostream>

int main() {
  Simulation3D s("map05.dat", 1, 0.03, 5, 11);
  
  s.saveV("0.dat");
  std::cout << "Done 0/4\n";

  for (int i = 0; i < 2000; ++i) 
    s.evolve();
  s.saveV("20.dat");
  std::cout << "Done 1/4\n";

  for (int i = 0; i < 12000; ++i)
    s.evolve();
  s.saveV("140.dat");
  std::cout << "Done 2/4\n";

  for (int i = 0; i < 86000; ++i)
    s.evolve();
  s.saveV("1000.dat");
  std::cout << "Done 3/4\n";
  
  for (int i = 0; i < 100000; ++i)
    s.evolve();
  s.saveV("2000.dat");
  std::cout << "Done 4/4\n";

  std::cout << "Fine\n";
  return 0;
}