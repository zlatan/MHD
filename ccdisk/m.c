#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "parameters.h"

void set_initial_values_linear_term(LinearTerm *lt,double DQx,double DQy,double DQz,double DVQ);

int main(void)
{
    double DQx=Qx0/Nx, DQy=Qy0/Ny, DQz=Qz0/Nz;
    double DVQ=DQx*DQy*DQz/((2*M_PI)*(2*M_PI)*(2*M_PI));

	NonlinearTerm bb,vv,vb,bv;
    LinearTerm v,b,fv,fb;

    set_initial_values_linear_term(&v,DQx,DQy,DQz,DVQ);
    set_initial_values_linear_term(&b,DQx,DQy,DQz,DVQ);
}