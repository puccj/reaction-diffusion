//Discretization of second spatial derivatives
double secDerU(double h, int i, int j, int k) {
  return ( u[i+1][j][k] + u[i-1][j][k] + u[i][j+1][k] + u[i][j-1][k] + u[i][j][k+1] + u[i][j][k-1] - 6*u[i][j][k] ) / (h*h);
}
double secDerV(double h, int i, int j, int k) {
  return ( v[i+1][j][k] + v[i-1][j][k] + v[i][j+1][k] + v[i][j-1][k] + v[i][j][k+1] + v[i][j][k-1] - 6*v[i][j][k] ) / (h*h);
}

  //domain intervals of the 3D domain embedding the narrow band domain
  double a,b;
  double c,d;
  double e,f;

  int Nx, Ny, Nz;
  //uniform mesh size
  double h = (b-a) / Nx; // = (d-c)/Ny = (f-e)/Nz

  /*
  //discrete domain (3D array)
  double domain[(int)h*Nx][(int)h*Ny][(int)h*Nz];

  //approximation of u(xi, yj, zk, n/\t) and v(xi, yj, zk, n/\t)
  double u[(int)h*Nx][(int)h*Ny][(int)h*Nz];
  double v[(int)h*Nx][(int)h*Ny][(int)h*Nz];
  */
  
  
  double u[10][10][10];
  double v[10][10][10];

int main() {
  //Final time
  double T;
  //Total number of time steps
  int Nt;
  //Time step 
  double DeltaT = T/Nt;

  //Fast and accurate numerical method for motion by mean curvature of curves on a surfice in 3D space
  //Discretization of the reaction-diffusion system using explicit scheme:

  // (6) (7)

  //Discretization of second spatial derivative

  //Diffusion coefficients
  double Du = 1;
  double Dv = 0.02;

  // Constants related to feed concentration
  double k1 = 9;
  double k2 = 11;

  // Discretization of the reaction diffusion system (4) and (5)
  
  // Numerical closest point (8) 
  //cp(int i, int, j, int k);

  // If point is not lying on a given computational grid, we obtain u^n(cp(i,j,k)) by computating trilinear interpolation
  // for fast computation, tabulated the interpolation stencil and three fractions for each boundary point

  /*
  //Omega = [0,10]x[0,10]x[0,10]
  int a = 0;
  int b = 10;
  int c = 0;
  int d = 10;

  int Du = 1;
  int k2 = 11;

  //mesh grid 101 x 101
  double grid[101][101];
  double h = 0.1;
  double DeltaT = 0.1* h*h;

  //Periodic bounding conditions

  //Initial condition
  double uBar = 1 + 0.04*k2*k2;
  double vBar = 0.2*k2;

  std::random_device seed;
  std::uniform_real_distribution<double> dist(-1,1);
  double 
  */

}