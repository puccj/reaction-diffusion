//#include "matrix.h"
#include "simulation.h"
#include <iostream>
#include <random>
#include <fstream>

int main () {
  /*
  double a = m[3][2];
  std::cout << a << '\n';

  a = -3;

  m[2][3] = a;

  m.print();

  */

  double h = 0.1;
  double dt = 0.1*h*h;
  Simulation s(0,10,0,10, h);
  s.setDu(1);
  s.setk2(11);
  s.setDv(0.02);
  s.setk1(9);
  
  int barWidth = 70;
  int evolutions = 1000;
  std::cout << "Evolving...\n";
  for (int i = 0; i < evolutions; ++i) {
    s.evolve(dt);

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
  }

  //delete progress bar
  std::cout << "Done!";
  for (int i = 0; i < barWidth +2; ++i)
    std::cout << ' ';
  std::cout << '\n';

  s.saveV("simulationA.dat");


  return 0;
}