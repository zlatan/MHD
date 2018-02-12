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


	double rnvx[2*Nx+2][2*Ny+2][2*Nz+2], rnvy[2*Nx+2][2*Ny+2][2*Nz+2], rnvz[2*Nx+2][2*Ny+2][2*Nz+2];
	double rnbx[2*Nx+2][2*Ny+2][2*Nz+2], rnby[2*Nx+2][2*Ny+2][2*Nz+2], rnbz[2*Nx+2][2*Ny+2][2*Nz+2];
	int Lx=2*Nx, Ly= 2*Ny, Lz=2*Nz;
	double rnx, rny, rnz, vn, bn, Q;

	for(int ix=0;ix<=Lx; ix++)
	{double Qx=-Qx0+ix*DQx;
		for(int iy=0;iy<=Ly; iy++)
		{double Qy=-Qy0+iy*DQy;
			for(int iz=0;iz<=Lz; iz++)
			{double Qz=-Qz0+iz*DQz;
      
				rnvx[ix][iy][iz]=(vv.xx[ix][iy][iz]+bb.xx[ix][iy][iz])*Qx+(vv.xy[ix][iy][iz]+bb.xy[ix][iy][iz])*Qy+(vv.xz[ix][iy][iz]+bb.xz[ix][iy][iz])*Qz;
				rnvy[ix][iy][iz]=(vv.yx[ix][iy][iz]+bb.yx[ix][iy][iz])*Qx+(vv.yy[ix][iy][iz]+bb.yy[ix][iy][iz])*Qy+(vv.yz[ix][iy][iz]+bb.yz[ix][iy][iz])*Qz;
				rnvz[ix][iy][iz]=(vv.zx[ix][iy][iz]+bb.zx[ix][iy][iz])*Qx+(vv.zy[ix][iy][iz]+bb.zy[ix][iy][iz])*Qy+(vv.zz[ix][iy][iz]+bb.zz[ix][iy][iz])*Qz;

				rnbx[ix][iy][iz]=(bv.xx[ix][iy][iz]-vb.xx[ix][iy][iz])*Qx+(bv.xy[ix][iy][iz]-vb.xy[ix][iy][iz])*Qy+(bv.xz[ix][iy][iz]-vb.xz[ix][iy][iz])*Qz;
				rnby[ix][iy][iz]=(bv.yx[ix][iy][iz]-vb.yx[ix][iy][iz])*Qx+(bv.yy[ix][iy][iz]-vb.yy[ix][iy][iz])*Qy+(bv.yz[ix][iy][iz]-vb.yz[ix][iy][iz])*Qz;
				rnbz[ix][iy][iz]=(bv.zx[ix][iy][iz]-vb.zx[ix][iy][iz])*Qx+(bv.zy[ix][iy][iz]-vb.zy[ix][iy][iz])*Qy+(bv.zz[ix][iy][iz]-vb.zz[ix][iy][iz])*Qz;

		         	Q=sqrt(Qx*Qx+Qy*Qy+Qz*Qz); if (Q==0.0) continue;
				rnx=Qx/Q; rny= Qy/Q; rnz= Qz/Q;
				vn=rnvx[ix][iy][iz]*rnx+ rnvy[ix][iy][iz]*rny+rnvz[ix][iy][iz]*rnz;
				bn=rnbx[ix][iy][iz]*rnx+ rnby[ix][iy][iz]*rny+rnbz[ix][iy][iz]*rnz;
				rnvx[ix][iy][iz]-=rnx*vn; rnvy[ix][iy][iz]-=rny*vn; rnvz[ix][iy][iz]-=rnz*vn;
				rnbx[ix][iy][iz]-=rnx*bn; rnby[ix][iy][iz]-=rny*bn; rnbz[ix][iy][iz]-=rnz*bn;
				printf("%g %g %g\n",Qy, Qz, vv.xx[ix][iy][iz]);
			}
		}
	}



   return 0;
}

