#include "variables.h"
#include "parameters.h"
#include <math.h>

void set_initial_values_linear_term(LinearTerm *lt,  double DQx,double DQy,double DQz,double DVQ)
{
    for (int ix=0;ix<=2*Nx; ix++)
	    for (int iy=0;iy<=2*Ny; iy++)
	        for (int iz=0;iz<=2*Nz; iz++)
                {
                double Qx=-Qx0+ix*DQx;
                double Qy=-Qy0+iy*DQy;
                double Q2 = Qx * Qx + Qy * Qy;
                double norm = 1 / sqrt(Q2);

                if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny)) == 0) 
                    continue;
                

                lt->x[ix][iy][iz] = vphi * -Qy * norm;
                lt->y[ix][iy][iz] = vphi * Qx * norm;
                lt->z[ix][iy][iz] = 0;
                }

    lt->x[Nx + 1][Ny][Nz] = 0.0;
    lt->y[Nx + 1][Ny][Nz] = vc;

    lt->x[Nx][Ny + 1][Nz] = -vc;
    lt->y[Nx][Ny + 1][Nz] = 0.0;

    lt->x[Nx - 1][Ny][Nz] = 0.0;
    lt->y[Nx - 1][Ny][Nz] = -vc;

    lt->x[Nx][Ny - 1][Nz] = vc;
    lt->y[Nx][Ny - 1][Nz] = 0.0;
}