#include <stdio.h>
#include <math.h>

#define D 200

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

Function f;
int i=0;

double a,b;
int N,Ni;
double xi,dx;

Ni=6;
a=-3.14;
b=3.14;
N=100;
dx=(b-a)/N;

/*
double dg=M_PI/180.;
int i, N=4;
function.argument[0]=0.;
function.argument[1]=30.;
function.argument[2]=45.;
function.argument[3]=60.;
function.argument[4]=90.;

for(i=0;i<=4;i++) 
  function.value[i]=sin(function.argument[i]*dg);
*/

//printf("f[%i]=%g\n",N,aitken(&function,N,50.0));


//for(double y=-3.14;y<=3.14;y=y+0.1)
//{



for(i=0;i<=N;i++)
{
   xi = a + (dx*i);
   f.argument[i] = xi;
   f.value[i] = exp(-xi*xi);
}

double x[Ni], fi[Ni];

//printf("%g %g\n",x,exp(-x*x));
 
//printf("%g %g\n",y,aitken_pointer(&function,62,y));
//}

double xx=b;

double ff;
int j;


j=(int) ( (xx-a)/dx );
printf("j(0)=%i\n",j);
j-=2;
if (j<1) 
    j=1;

if (j + Ni > N) 
   j = N - Ni;

printf("j=%i N=%i Ni=%i\n",j,N,Ni);

for(i=0;i<=Ni;i++)
{
x[i]=f.argument[i+j];
fi[i]=f.value[i+j];
printf("x[%i]=%g, fi[%i]=%g\n",i,x[i],i,fi[i]);
}// prezarezhdane sys kys masiv


for(int k=0;k<Ni; k++)
 {
 for(i=k+1;i<=Ni;i++) 
  {
   fi[i]=( (x[k]-xx)*fi[i] - (x[i]-xx)*fi[k] )/( (x[k]-xx)- (x[i]-xx) );
  }
ff=fi[k+1];
 }
printf("ff=%g, xx=%g\n",ff,xx);

// printf("%g %g\n",0.1234,aitken(&f,N,0.1234));


//for(double y=0;y<=3.14;y=y+0.1){
//}




return 0;
}
