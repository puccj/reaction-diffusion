//#include "matrix.h"
#include "simulation.h"
#include <iostream>
#include <random>
#include <fstream>

int main () {
  double h = 0.1;
  double dt = 0.1*h*h;
  Simulation s(0,10,0,10, h);
  s.setDu(1);
  s.setk2(11);
  s.setDv(0.02);
  s.setk1(9);
  
  int barWidth = 70;
  int evolutions = 14000;
  std::cout << "Evolving...\n";
  s.evolve(dt, false);
  for (int i = 0; i < 2000; ++i) {
    s.evolve(dt, true);

    /*
    //print progress bar
    std::cout << "[";
    int pos = barWidth * i/evolutions;
    for (int j = 0; j < barWidth; ++j) {
      if (j < pos) std::cout << '=';
      else if (j == pos) std::cout << '>';
      else std::cout << ' ';
    }
    std::cout << "] " << int(i* 100/evolutions) << " %\r";
    std::cout.flush();
    */
  }
  std::cout << "Step 1 done.\n";
  s.evolve(dt, false);

  for (int i = 0; i < 12000; ++i) {
    s.evolve(dt, true);
  }
  std::cout << "Step 2 done.\n";
  s.evolve(dt, false);

  for (int i = 0; i < 86000; ++i) {
    s.evolve(dt, true);
  }
  std::cout << "Step 3 done.\n";
  s.evolve(dt, false);

  for (int i = 0; i < 100000; ++i) {
    s.evolve(dt, true);
  }
  std::cout << "Step 4 done.\n";

  //delete progress bar
  std::cout << "Done!";
  /*
  for (int i = 0; i < barWidth +2; ++i)
    std::cout << ' ';
  std::cout << '\n';
  */

  s.saveV("stripes.dat");

  return 0;
}