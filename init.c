#include "parameters.h"
#include "variables.h"
#include <stdlib.h>
#include <math.h>

void set_initial_values_nonlinear_term(NonlinearTerm *nt)
{
    for (int ix=0;ix<=2*Nx; ix++)
	for (int iy=0;iy<=2*Ny; iy++)
	    for (int iz=0;iz<=2*Nz; iz++)
	    {
	      nt->xx[ix][iy][iz]=(rand()-.5)*1E-10;
	      nt->xy[ix][iy][iz]=(rand()-.5)*1E-10;
	      nt->xz[ix][iy][iz]=(rand()-.5)*1E-10;
	      nt->yx[ix][iy][iz]=(rand()-.5)*1E-10;
	      nt->yy[ix][iy][iz]=(rand()-.5)*1E-10;
	      nt->yz[ix][iy][iz]=(rand()-.5)*1E-10;
      	      nt->zx[ix][iy][iz]=(rand()-.5)*1E-10;
      	      nt->zy[ix][iy][iz]=(rand()-.5)*1E-10;
      	      nt->zz[ix][iy][iz]=(rand()-.5)*1E-10;
	    }
}

void set_initial_values_linear_term(LinearTerm *lt,  double DQx,double DQy,double DQz,double DVQ)
{
    for (int ix=0;ix<=2*Nx; ix++)
	{double Qx=-Qx0+ix*DQx;
	for (int iy=0;iy<=2*Ny; iy++)
	    {double Qy=-Qy0+iy*DQy;
	    for (int iz=0;iz<=2*Nz; iz++)
		{double Qz=-Qz0+iz*DQz;
		 double Q2;
                 double Q;
		 double rnx, rny, rnz, n;

                 lt->x[ix][iy][iz]=1.0;
		 lt->y[ix][iy][iz]=1.0;
		 lt->z[ix][iy][iz]=0.0;

 		 Q2=Qx*Qx+Qy*Qy+Qz*Qz;
   		 Q=sqrt(Q2);

		if ( ((ix - Nx)*(ix - Nx) + (iy-Ny)*(iy-Ny) + (iz-Nz)*(iz-Nz)) == 0 )
			continue;

		rnx=Qx/Q; rny=Qy/Q; rnz=Qz/Q;
                n = lt->x[ix][iy][iz]*rnx + lt->y[ix][iy][iz]*rny + lt->z[ix][iy][iz]*rnz;
		lt->x[ix][iy][iz]-= n*rnx;
		lt->y[ix][iy][iz]-= n*rny;
		lt->z[ix][iy][iz]-= n*rnz;
	    }
	}
     }
}

void set_DV(double DVQ)
{
      for (int ix=0;ix<=2*Nx; ix++)
	  {
	  double Qx=-Qx0+ix*DQx;
	for (int iy=0;iy<=2*Ny; iy++)
	    {
	     double Qy=-Qy0+iy*DQy;
	    for (int iz=0;iz<=2*Nz; iz++)
     		{
		 double Qz=-Qz0+iz*DQz;
	         DV[ix][iy][iz]=DVQ;
		}
	    }
	  }

}


void set_initial_values_linear_term_gauss(LinearTerm *lt,double DQx,double DQy,double DQz,double DVQ)
{
  printf("Nx=%g\n",Nx);
    for(int ix=0;ix<=2*Nx; ix++)
	{
	double Qx=-Qx0+ix*DQx;
	for(int iy=0;iy<=2*Ny; iy++)
	    {
            double Qy=-Qy0+iy*DQy;
	    for(int iz=0;iz<=2*Nz; iz++)
	       {
		 double Qz=-Qz0+iz*DQz;
		 double Q2;
                      Q2=Qx*Qx + Qy*Qy + Qz*Qz;
		      lt->x[ix][iy][iz]=exp(-Q2);
		      lt->y[ix][iy][iz]=1.0;
		      lt->z[ix][iy][iz]=0.0;

	}
     }
   }
}


