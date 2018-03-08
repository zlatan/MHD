#include <stdio.h>
#include <math.h>

#define D 4

typedef struct Function
{
      double argument[D+1];
      double value[D+1];
} Function;



double aitken(Function *function, int N, double point)
{
// TODO: SORT
for(int k=0;k<N; k++)
 for(int i=k+1;i<=N;i++) 
  function->value[i]=((function->argument[k]-point)*function->value[i] - (function->argument[i]-point)*function->value[k])
			 / ((function->argument[k]-point) - (function->argument[i]-point));

return function->value[N];
}


int main(void)
{

Function function;
double dg=M_PI/180.;
int i, N=4;

function.argument[0]=0.;
function.argument[1]=30.;
function.argument[2]=45.;
function.argument[3]=60.;
function.argument[4]=90.;


for(i=0;i<=4;i++) 
  function.value[i]=sin(function.argument[i]*dg);


printf("f[%i]=%g\n",N,aitken(&function,N,50.0));

return 0;
}
