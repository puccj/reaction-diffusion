#ifndef SIMULATION3D_H
#define SIMULATION3D_H

#include "surface.h"
//#include "utils.hpp"

class Simulation3D {
 private:
  //map a surface point to its 6 neighbors + itself.
  //Each of the neighbors 'n' is in turn mapped to:
  // - pair.fisrt:  the indexes (inside of _s array) of the 8 points surrounding closest(n)
  // - pair.second: the weights of the interpolation with each of these 8 points
  struct Map {
    std::pair<int[8], double[8]> neighbor[7]; 
  };

  int _nPoints;
  Surface _s;
  double* _u;   //concentration of the inhibitor
  double* _v;   //concentration of the activator
  Map* _map;    //map each point to the interpolations of the closests point of its neighbour
  //u, v and map have lenght nPoints

  double _h;      //Uniform mesh size
  double _delta;  //half-thickness of narrow band domain
  double _Du;
  double _Dv;
  double _k1;
  double _k2;

 public:
  Simulation3D(Interval x = {0,10}, Interval y = {0,10}, Interval z = {0,10}, double h = 0.1, double k2 = 11);
  Simulation3D(std::string mapFile, double Du, double Dv, double k1, double k2 = 11);
  //Simulation3D(double iMin = 0, double iMax = 10, double jMin = 0, double jMax = 10, double kMin, double kMax, double h = 1, double k2 = 11);
  ~Simulation3D();

  void setDu(double val) {_Du = val; };
  void setDv(double val) {_Dv = val; };
  void setk1(double val) {_k1 = val; };
  void setk2(double val);

  //save the map in a file. Format: 7 rows for each point (one per neighbor 'n' + itself).
  //Each row has 8 + 8 value for indexes and weights of the 8 point surrounding closest(n)
  void saveMap(std::string filename = "map.dat");

  //save surface points in a file. Format: 3 cols of data: x y z
  void saveSurface(std::string filename = "surface.dat");

  //save activator concentration in a file. Format: 1 col of data
  void saveV(std::string filename = "value.dat");

  //calculate the next step of the Tensor, given the time interval dt
  void evolve(double dt = 0);

 private:
  //return phi (signed distance function) calculated at x,y,z
  double distFunc(double x, double y, double z);

  //create the first UV value
  void createUV();

  //set the first values of _u and _v
  void fillFirstUV();

  //create the map between each surface point 'p' and the interpolation(closest(neighbour(p)))
  void constructMap();

  //return the point on the surface which is closest to given point
  Point closest(Point const& p);
  
  //load the map between each surface point 'p' and the interpolation(closest(neighbour(p))).
  //Returns false if file is not found.
  bool loadMap(std::string filename = "map.dat");

  //return the approximation of the second derivative of given concentration at point corresponding to index
  double der2(char concentration, int pointIndex);  
};

#endif //SIMULATION3D_H