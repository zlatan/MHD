#include "parameters.h"
#include "variables.h"
#include <stdlib.h>

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

