#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "parameters.h"
#include "main.h"

int main(void)
{
        double DQx=Qx0/Nx, DQy=Qy0/Ny, DQz=Qz0/Nz;
	double DVQ=DQx*DQy*DQz/((2*M_PI)*(2*M_PI)*(2*M_PI));
	double Tempi=10.; 
        double Nt=2,it=0,t=0;
        double dt=Tempi/Nt;

	double Ei[2*Nx][2*Nx][2*Nx];
	double Ef[2*Nx][2*Nx][2*Nx];

	NonlinearTerm bb,vv,vb,bv;
	LinearTerm v,b,fv,fb;

	set_initial_values_nonlinear_term_zero(&bb);
	set_initial_values_nonlinear_term_zero(&vv);
	set_initial_values_nonlinear_term_zero(&vb);
	set_initial_values_nonlinear_term_zero(&bv);

        set_initial_values_linear_term(&v,DQx,DQy,DQz,DVQ);
        set_initial_values_linear_term(&b,DQx,DQy,DQz,DVQ);


        calculate_nonlinear_term(&bb, &vv, &bv,  &vb, &b, &v, DQx, DQy, DQz, DVQ);
	calculate_linear_term(&bb, &vv, &bv,  &vb, &b, &v, &fv, &fb, DQx, DQy, DQz, DVQ);
	integrate_over_t(&b, &v, &fv, &fb, dt);
	
	for(int ix=0;ix<=2*Nx; ix++)
		for(int iy=0;iy<=2*Ny; iy++)
			for(int iz=0;iz<=2*Nz; iz++)
			Ei[ix][iy][iz] = v.x[ix][iy][iz]*v.x[ix][iy][iz] 
					+ v.y[ix][iy][iz]*v.y[ix][iy][iz] 
					+ v.z[ix][iy][iz]*v.z[ix][iy][iz]
					+ b.x[ix][iy][iz]*b.x[ix][iy][iz] 
					+ b.y[ix][iy][iz]*b.y[ix][iy][iz] 
					+ b.z[ix][iy][iz]*b.z[ix][iy][iz];



        calculate_nonlinear_term(&bb, &vv, &bv,  &vb, &b, &v, DQx, DQy, DQz, DVQ);
	calculate_linear_term(&bb, &vv, &bv,  &vb, &b, &v, &fv, &fb, DQx, DQy, DQz, DVQ);
	integrate_over_t(&b, &v, &fv, &fb, dt);

	for(int ix=0;ix<=2*Nx; ix++)
		for(int iy=0;iy<=2*Ny; iy++)
			for(int iz=0;iz<=2*Nz; iz++)
			Ef[ix][iy][iz] = v.x[ix][iy][iz]*v.x[ix][iy][iz] 
					+ v.y[ix][iy][iz]*v.y[ix][iy][iz] 
					+ v.z[ix][iy][iz]*v.z[ix][iy][iz] 
					+ b.x[ix][iy][iz]*b.x[ix][iy][iz] 
					+ b.y[ix][iy][iz]*b.y[ix][iy][iz] 
					+ b.z[ix][iy][iz]*b.z[ix][iy][iz];


	for(int ix=0;ix<=2*Nx; ix++)
	   {
 	   double Qx=-Qx0+ix*DQx;
			for(int iz=0;iz<=2*Nz; iz++)
			   {
			   double Qz=-Qz0+iz*DQz;
			   printf("%g %g %g\n",Qx,Qz,pow( (1/(2.0*Tempi))*log(Ef[ix][Ny][iz]/Ei[ix][Ny][iz]),2));
			}
	}



   return 0;
}

