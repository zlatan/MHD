#include <stdio.h>
#include <math.h>
#include "variables.h"
#include "parameters.h"
#include "init.h"

int main(void)
{
	
	NonlinearTerm nv;
	set_initial_values_nonlinear_term(&nv);
        printf( "value=%g\n", nv.xx[1][1][1]);

return 0;
}



