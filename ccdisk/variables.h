#include "parameters.h"

typedef struct NonlinearTerm 
{
      double xx[2*Nx+2][2*Ny+2][2*Nz+2];
      double xy[2*Nx+2][2*Ny+2][2*Nz+2];
      double xz[2*Nx+2][2*Ny+2][2*Nz+2];
      double yx[2*Nx+2][2*Ny+2][2*Nz+2];
      double yy[2*Nx+2][2*Ny+2][2*Nz+2];
      double yz[2*Nx+2][2*Ny+2][2*Nz+2];
      double zx[2*Nx+2][2*Ny+2][2*Nz+2];
      double zy[2*Nx+2][2*Ny+2][2*Nz+2];
      double zz[2*Nx+2][2*Ny+2][2*Nz+2];
} NonlinearTerm;

typedef struct LinearTerm
{
      double x[2*Nx+2][2*Ny+2][2*Nz+2];
      double y[2*Nx+2][2*Ny+2][2*Nz+2];
      double z[2*Nx+2][2*Ny+2][2*Nz+2];
} LinearTerm;

typedef struct vector
{
      double x;
      double y;
      double z;
} vector;

typedef struct unit_vector
{
      vector phy;
      vector r;
      vector theta;
} unit_vector;