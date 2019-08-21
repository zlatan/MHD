#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "parameters.h"

void set_initial_values_linear_term(LinearTerm *lt,double DQx,double DQy,double DQz,double DVQ);
void set_initial_values_nonlinear_term_zero(NonlinearTerm *nt);
void calculate_nonlinear_term(NonlinearTerm *bb, NonlinearTerm *vv, NonlinearTerm *bv, NonlinearTerm *vb,LinearTerm *b, LinearTerm *v, double DQx,double DQy,double DQz,double DVQ);
void calculate_linear_term(NonlinearTerm *bb, NonlinearTerm *vv, NonlinearTerm *bv, NonlinearTerm *vb, LinearTerm *b, LinearTerm *v, double DQx,double DQy,double DQz,double DVQ);

int main(void)
{
    double DQx=Qx0/Nx, DQy=Qy0/Ny, DQz=Qz0/Nz;
    double DVQ=DQx*DQy*DQz/((2*M_PI)*(2*M_PI)*(2*M_PI));
    double E = 0.0;

	NonlinearTerm bb,vv,vb,bv;
    LinearTerm v,b;

    set_initial_values_linear_term(&v,DQx,DQy,DQz,DVQ);
    set_initial_values_linear_term(&b,DQx,DQy,DQz,DVQ);

    // set_initial_values_nonlinear_term_zero(&bb);
	// set_initial_values_nonlinear_term_zero(&vv);
	// set_initial_values_nonlinear_term_zero(&vb);
	// set_initial_values_nonlinear_term_zero(&bv);


for (int i=0; i< 1000 ; i ++){
    set_initial_values_nonlinear_term_zero(&bb);
	set_initial_values_nonlinear_term_zero(&vv);
	set_initial_values_nonlinear_term_zero(&vb);
	set_initial_values_nonlinear_term_zero(&bv);

    calculate_nonlinear_term(&bb, &vv, &bv, &vb, &b, &v, DQx, DQy, DQz, DVQ);
       calculate_linear_term(&bb, &vv, &bv, &vb, &b, &v, DQx, DQy, DQz, DVQ);

    for (int ix = 0; ix <= 2*Nx; ix++) {
        for (int iy = 0; iy <= 2*Ny; iy++) {
            for (int iz = 0; iz <= 2*Nz; iz++) {
                 E += (v.x[ix][iy][iz] * v.x[ix][iy][iz] + v.y[ix][iy][iz] * v.y[ix][iy][iz] + v.z[ix][iy][iz] * v.z[ix][iy][iz] 
                    +  b.x[ix][iy][iz] * b.x[ix][iy][iz] + b.y[ix][iy][iz] * b.y[ix][iy][iz] + b.z[ix][iy][iz] * b.z[ix][iy][iz] )*DVQ;
            }
        }
    }

  }

}