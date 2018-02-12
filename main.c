#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "parameters.h"
#include "init.h"
#include "nonlinear.h"

int main(void)
{
        double DQx=Qx0/Nx, DQy=Qy0/Ny, DQz=Qz0/Nz;
	double DVQ=DQx*DQy*DQz/((2*M_PI)*(2*M_PI)*(2*M_PI));

	NonlinearTerm bb,vv,vb,bv;
	LinearTerm v,b,fv,fb;

	set_initial_values_nonlinear_term_zero(&bb);
	set_initial_values_nonlinear_term_zero(&vv);
	set_initial_values_nonlinear_term_zero(&vb);
	set_initial_values_nonlinear_term_zero(&bv);
        
	set_initial_values_linear_term_gauss(&v,DQx,DQy,DQz,DVQ);
        calculate_nonlinear_term(&bb, &vv, &bv,  &vb, &b, &v, DQx, DQy, DQz, DVQ);


	for(int ix=0;ix<=2*Nx; ix++)
	{double Qx=-Qx0+ix*DQx;
		for(int iy=0;iy<=2*Ny; iy++)
		{double Qy=-Qy0+iy*DQy;
			for(int iz=0;iz<=2*Nz; iz++)
			{double Qz=-Qz0+iz*DQz;
			printf("%g %g %g %g\n",Qx, Qy, Qz, vv.xx[ix][iy][iz]);
				//printf("%g %g %g\n",Qy, Qz, vv.xx[ix][iy][iz]);
			}
		}
	}



   return 0;
}

