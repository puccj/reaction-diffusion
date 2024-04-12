#include "simulation3.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>

Simulation3D::Simulation3D(Interval x, Interval y, Interval z, double h, double k2) 
  : _nPoints{0},
    _u{nullptr},
    _v{nullptr},
    _map{nullptr},
    _h{h},
    _delta{1.1 * sqrt(3) * h},
    _k2{k2} 
{
  int domainPoints = (x.max - x.min)*(y.max - y.min)*(z.max - z.min) / (h*h*h);
  std::cout << "Debug: domainPoints = " << domainPoints << '\n';

  // Create surface and get _nPoints
  Point* temp = new Point[domainPoints];

  //for every point of the descrete domain..
  for (double i = x.min; i < x.max; i+=h) {
    for (double j = y.min; j < y.max; j+=h) {
      for (double k = z.min; k < z.max; k+=h) {
        //..add it to the surface definition if it's inside the narrow-band domain
        double dist = phi(i,j,k);
        if (dist > -_delta && dist < _delta) {
          temp[_nPoints] = {i,j,k};
          ++_nPoints;
        }
      }
    }
  }
  
  std::cout << "Debug: surfacePoints = " << _nPoints << '\n';
  _s = Surface(_nPoints, temp);  //call copy constructor, so temp is safe to be deleted
  delete[] temp;

  // Create uv and map
  createUV();
  constructMap();
}

Simulation3D::Simulation3D(std::string mapFile, double Du, double Dv, double k1, double k2) 
  : _s{Surface()},
    _u{nullptr},
    _v{nullptr},
    _map{nullptr},
    _Du{Du},
    _Dv{Dv},
    _k1{k1},
    _k2{k2}
{
  if (!loadMap(mapFile))
    throw std::runtime_error{"loadMap() cannot open file."};
  createUV();
}

Simulation3D::~Simulation3D() {
  delete[] _u;
  delete[] _v;
  delete[] _map;
}

void Simulation3D::setk2(double val) {
  _k2 = val;
  //re-create the first matrixes, since they depend on k2 value
  fillFirstUV();
}

void Simulation3D::saveMap(std::string filename) {
  //format: 7 rows for each point (one row for each neighbor 'n' + itself)
  //        in each of these rows: index and weight of the 8 points surrounding n

  std::cout << "Saving map...     ";
  std::fstream fout(filename, std::ios::out);

  fout << _h << ' ' << _nPoints << "\n\n";
  
  for (int i = 0; i < _nPoints; ++i) {
    for (int n = 0; n < 7; ++n) {
      for (int j = 0; j < 8; ++j) {
        fout << _map[i].neighbor[n].first[j] << ' ';
        fout << _map[i].neighbor[n].second[j] << ' ';
      }
      fout << '\n';
    }
    fout << "\n";
  }

  std::cout << "Done!\n";
}

void Simulation3D::saveSurface(std::string filename) {
  std::cout << "Saving surface...     ";
  if (_s.nPoints() == 0) {
    std::cerr << "Error, cannot save surface of a simulation created from a map.\n";
    return;
  }

  std::fstream fout(filename, std::ios::out);
  fout << _s;
  fout.close();

  std::cout << "Done!\n";
}

void Simulation3D::saveU(std::string filename) {
  std::fstream fout(filename, std::ios::out);
  for (int i = 0; i < _nPoints; ++i)
    fout << _u[i] << '\n';

  fout.close();
}

void Simulation3D::saveV(std::string filename) {
  //format: 1 col of data: value of V at each point
  std::fstream fout(filename, std::ios::out);
  for (int i = 0; i < _nPoints; ++i)
    fout << _v[i] << '\n';

  fout.close();
}

void Simulation3D::evolve(double dt) {
  if (dt == 0)
    dt = 0.1*_h*_h*_h;
  
  double* nextU = new double[_nPoints];
  double* nextV = new double[_nPoints];

  for (int i = 0; i < _nPoints; ++i) {
    //u_ijk is not _u[i], because you need to take the interpolation of the closest point of i
    double u_ijk = 0;
    double v_ijk = 0;

    std::pair<int[8], double[8]> stencil = _map[i].neighbor[0];
    for (int j = 0; j < 8; ++j) {
      int index = stencil.first[j];
      double weight = stencil.second[j];
      u_ijk += _u[index] * weight;
      v_ijk += _v[index] * weight;
    }

    //calculate evolution
    double frac = u_ijk*v_ijk/ (1+ v_ijk*v_ijk);
    nextU[i] = u_ijk + (_Du*der2('u',i) + _k1 * frac) *dt;
    nextV[i] = v_ijk + (_Dv*der2('v',i) + _k2 - v_ijk  - 4*frac) *dt;
  }

  //release data from u and v
  delete[] _u;
  delete[] _v;

  //move data from next
  _u = nextU;
  _v = nextV;

  //just for safety
  nextU = nullptr;
  nextV = nullptr;
}


//private


double Simulation3D::phi(double x, double y, double z) {
  //default function: a sphere centred in (5,5,5) with r=4.5
  return std::sqrt((x-5)*(x-5) + (y-5)*(y-5) + (z-5)*(z-5)) - 4.5;
}

void Simulation3D::createUV() {
  if (_u != nullptr) {
    std::cerr << "Warning: u already created. Deleting it before recreating\n";
    delete[] _u;
  }
  if (_v != nullptr) {
    std::cerr << "Warning: v already created. Deleting it before recreating\n";
    delete[] _v;
  }

  _u = new double[_nPoints];
  _v = new double[_nPoints];

  fillFirstUV();
}

void Simulation3D::fillFirstUV() {
  std::random_device seed;
  std::cout << "Debug: seed = " << seed() << '\n';
  std::mt19937 gen(seed());
  std::uniform_real_distribution<double> dist(-0.1,0.1);

  for (int i = 0; i < _nPoints; ++i) {
    //default functions
    _u[i] = 1+ 0.04*_k2*_k2 + dist(gen);
    _v[i] = 0.2*_k2 + dist(gen);
  }
}

void Simulation3D::constructMap() {
  auto start_time = std::chrono::high_resolution_clock::now();
  std::cout << "Constructing map...     ";
  _map = new Map[_nPoints];

  //std::fstream fout("debug.txt", std::ios::out);

  //for each surface point...
  for (int i = 0; i < _nPoints; ++i) {
    //fout << "Debug: Point " << i << '\n';

    //get the coordinate of the surface point
    Point p = _s[i];
    
    //fout << "Debug: Coordinates: " << p << '\n';

    //construct the 6 neighbor points + itself
    Point n[7];
    n[0] = {p.x,    p.y,    p.z   };
    n[1] = {p.x-_h, p.y,    p.z   };
    n[2] = {p.x+_h, p.y,    p.z   };
    n[3] = {p.x,    p.y-_h, p.z   };
    n[4] = {p.x,    p.y+_h, p.z   };
    n[5] = {p.x,    p.y,    p.z-_h};
    n[6] = {p.x,    p.y,    p.z+_h};
    
    /*for (int j = 0; j < 7; ++j)
      fout << "Debug: Neighbors: " << n[j] << '\n';
    */

    //for each of these neighbors..
    for (int j = 0; j < 7; ++j ) {
      
      //fout << "Debug: Neighbor " << j << ' ' << n[j] << ":\n";

      //find the point on the surface which is closest
      Point c = closest(n[j]);
      
      //fout << "Debug: Closest point: " << c << '\n';
      
      //(A)  find the indexes of the 8 surface points surrounding the closest point of neighbor[j]:
      //(A1) find the points (actually, only the up-left-front)
      double x0 = (int)(c.x/_h) *_h;
      double y0 = (int)(c.y/_h) *_h;
      double z0 = (int)(c.z/_h) *_h;

      //fout << "Debug: up-left-front point: " << Point{x0,y0,z0} << '\n';

      //(A2) find their indexes in _s
      for (int k = 0; k < _nPoints; ++k) {
        Point s = _s[k];
        
             if (s == Point{x0,    y0,    z0   })   _map[i].neighbor[j].first[0] = k;
        else if (s == Point{x0,    y0,    z0+_h})   _map[i].neighbor[j].first[1] = k;
        else if (s == Point{x0,    y0+_h, z0   })   _map[i].neighbor[j].first[2] = k;
        else if (s == Point{x0,    y0+_h, z0+_h})   _map[i].neighbor[j].first[3] = k;
        else if (s == Point{x0+_h, y0,    z0   })   _map[i].neighbor[j].first[4] = k;
        else if (s == Point{x0+_h, y0,    z0+_h})   _map[i].neighbor[j].first[5] = k;
        else if (s == Point{x0+_h, y0+_h, z0   })   _map[i].neighbor[j].first[6] = k;
        else if (s == Point{x0+_h, y0+_h, z0+_h})   _map[i].neighbor[j].first[7] = k;
      }

      /*for (int k = 0; k < 8; ++k)
        fout << "Debug: index of 8 points: " << _map[i].neighbor[j].first[k] << '\n';
      */

      //(B) find the weights of the interpolation
      double xd = (c.x - x0) /_h;
      double yd = (c.y - y0) /_h;
      double zd = (c.z - z0) /_h;

      _map[i].neighbor[j].second[0] = (1-xd) * (1-yd) * (1-zd);
      _map[i].neighbor[j].second[1] = (1-xd) * (1-yd) *   zd;
      _map[i].neighbor[j].second[2] = (1-xd) *   yd   * (1-zd);
      _map[i].neighbor[j].second[3] = (1-xd) *   yd   *   zd;
      _map[i].neighbor[j].second[4] =   xd   * (1-yd) * (1-zd);
      _map[i].neighbor[j].second[5] =   xd   * (1-yd) *   zd;
      _map[i].neighbor[j].second[6] =   xd   *   yd   * (1-zd);
      _map[i].neighbor[j].second[7] =   xd   *   yd   *   zd;

      /*for (int k = 0; k < 8; ++k)
        fout << "Debug: weight of 8 points: " << _map[i].neighbor[j].second[k] << '\n';
      */
    }
  }

  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time);
  std::cout << "Done! (" << duration.count() << " minutes)\n";
}

Point Simulation3D::closest(Point const& p) {
  double x = p.x;
  double y = p.y;
  double z = p.z;

  double gradient[3];   //nabla_h (phi)
  gradient[0] = (phi(x+_h, y   , z   ) - phi(x-_h, y,    z   )) /(2*_h);
  gradient[1] = (phi(x,    y+_h, z   ) - phi(x,    y-_h, z   )) /(2*_h);
  gradient[2] = (phi(x,    y,    z+_h) - phi(x,    y,    z-_h)) /(2*_h);

  //if gradient = 0, the point is already on the surface
  if (gradient[0] == 0 && gradient[1] == 0 && gradient[2] == 0) {
    return p;
  }

  // | nabla_h (phi) |
  double mod = std::sqrt(gradient[0]*gradient[0] + gradient[1]*gradient[1] + gradient[2]*gradient[2]);

  double dist_xyz = phi(x,y,z);

  Point r;
  r.x = x - gradient[0] * dist_xyz / mod;
  r.y = y - gradient[1] * dist_xyz / mod;
  r.z = z - gradient[2] * dist_xyz / mod;

  return r;
}

bool Simulation3D::loadMap(std::string filename) {
  std::cout << "Loading map...     ";
  std::fstream fin(filename, std::ios::in);
  if (!fin.is_open()) {
    std::cerr << "Error. Map file '"<< filename << "' couldn't be opened.\n";
    return false;
  }

  fin >> _h;
  fin >> _nPoints;

  _map = new Map[_nPoints];

  for (int i = 0; i < _nPoints; ++i) {
    for (int n = 0; n < 7; ++n) {
      for (int j = 0; j < 8; ++j) {
        fin >> _map[i].neighbor[n].first[j];
        fin >> _map[i].neighbor[n].second[j];
      }
    }
  }

  std::cout << "Done.\n";
  return true;
}

double Simulation3D::der2(char concentration, int pointIndex) {
  std::pair<int[8], double[8]> n[7] = {_map[pointIndex].neighbor[0],
                                       _map[pointIndex].neighbor[1],
                                       _map[pointIndex].neighbor[2],
                                       _map[pointIndex].neighbor[3],
                                       _map[pointIndex].neighbor[4],
                                       _map[pointIndex].neighbor[5],
                                       _map[pointIndex].neighbor[6]};

  //The function needs to return (down + up + right + left + front + back - 6*u_ijk) / (_h*_h);

  //interpolation: the value of function in a point p is a weighted sum of
  //the funcion value on the 8 points surrounding closest(p)

  //Remind that n[j].first  = indexes of the 8 points surrounding neighbor j
  //            n[j].second = weights of these 8 points in the interpolation

  double result = 0;

  if (concentration == 'u') {
    for (int i = 0; i < 8; ++i) {
      result -= 6 * _u[n[0].first[i]] * n[0].second[i];
      for (int j = 1; j < 7; ++j) {
        result += _u[n[j].first[i]] * n[j].second[i];
      }
    } 
  }
  else if (concentration == 'v') {
    for (int i = 0; i < 8; ++i) {
      result -= 6* _v[n[0].first[i]] * n[0].second[i];

      for (int j = 1; j < 7; ++j) {
        result += _v[n[j].first[i]] * n[j].second[i];
      } 
    }
  }
  else { 
    std::cerr << "Error: invalid concentration in der2 function\n";
    return nan("");
  }
  
  return result / (_h*_h);
}