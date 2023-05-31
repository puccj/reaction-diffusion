#include "simulation.h"
#include <iostream>

int main () {
  double h = 0.1;
  double dt = 0.1*h*h;

  for (double Dv = 0.01; Dv < 0.07; Dv += 0.01) {
    for (int k1 = 1; k1 < 12; ++k1) {

      Simulation s(0,10,0,10, h);
      s.setDu(1);
      s.setk2(11);
      s.setDv(Dv);
      s.setk1(k1);
      
      //for progress bar
      //int barWidth = 70;

      s.saveV("k1=" + std::to_string(k1) + ", dv=" + std::to_string(Dv) + " 0.dat");

      std::cout << "Evolving...\n";
      s.evolve(dt);
      for (int i = 0; i < 2000; ++i) {
        s.evolve(dt);

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
      s.saveV("k1=" + std::to_string(k1) + ", dv=" + std::to_string(Dv) + " 2000.dat");

      for (int i = 0; i < 12000; ++i)
        s.evolve(dt);

      std::cout << "Step 2 done.\n";
      s.saveV("k1=" + std::to_string(k1) + ", dv=" + std::to_string(Dv) + " 14000.dat");

      for (int i = 0; i < 86000; ++i)
        s.evolve(dt);
      std::cout << "Step 3 done.\n";
      s.saveV("k1=" + std::to_string(k1) + ", dv=" + std::to_string(Dv) + " 100000.dat");

      for (int i = 0; i < 100000; ++i)
        s.evolve(dt);
      std::cout << "Step 4 done.\n";
      s.saveV("k1=" + std::to_string(k1) + ", dv=" + std::to_string(Dv) + " 200000.dat");

      //delete progress bar
      std::cout << "Done!";
      /*
      for (int i = 0; i < barWidth +2; ++i)
        std::cout << ' ';
      std::cout << '\n';
      */
        
    }
  }

  std::cout << "Done All!\n";

  return 0;
}