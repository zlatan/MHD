#include <stdio.h>
#include <math.h>

#define D 10

typedef struct Function
{
      double argument[D+1];
      double value[D+1];
} Function;



double aitken(Function *f, int N, double point)
{
// TODO: SORT
for(int k=0;k<N; k++)
 for(int i=k+1;i<=N;i++) 
  f->value[i]=((f->argument[k]-point)*f->value[i] - (f->argument[i]-point)*f->value[k]) / ((f->argument[k]-point) - (f->argument[i]-point));
 //  f[i]      =                   ( (x[k]-xx)*f[i] - (x[i]-xx)*f[k] )                    /              ( (x[k]-xx)- (x[i]-xx) );
return f->value[N];
}

double calc(double x[],double f[],double xx, int Ni, int N, double a, double b)
{
Function function;
int i,j;
double dx;
double ff;
dx=(b-a)/N;
j=(int) ( (xx-a)/dx );

j-=2;
if (j<1) 
    j=1;

if (j+Ni > N) // Ni<N
   j = N - Ni;

for(i=0;i<=Ni;i++)
{
function.argument[i]=x[i+j];
function.value[i]=f[i+j];
}// lokalni interpolacinni tochki za stoinostta xx

ff=aitken(&function,Ni,xx);
return ff;
}


double shift(double x[],double f[],double xx, int Ni, int N, double a, double b, double s)
{
Function function;
int i,j;
double dx;
double ff;
dx=(b-a)/N;
j=(int) ( (xx-a)/dx );

j-=2;
if (j<1) 
    j=1;

if (j+Ni > N) // Ni<N
   j = N - Ni;

for(i=0;i<=Ni;i++)
{
function.argument[i]=x[i+j];
function.value[i]=f[i+j];
}// lokalni interpolacinni tochki za stoinostta xx

ff=aitken(&function,Ni,xx-s);
return ff;
}



int main(void)
{

double s=.3;
int i;
double a,b;
int N,Ni;
double xi,dx;

double xx,ff;
int j;
double result;
const double Pi=4.*atan(1.);
Ni=6;
a=-Pi;
b=Pi;
N=100;
dx=(b-a)/N;

double x[N+1],f[N+1], fn[N+1];

for(i=0;i<=N;i++)
{
   xi = a + (dx*i);
   x[i] = xi;
   f[i] = exp(-xi*xi);
}// tablichni stoinosti

//xx=0.1234;

for(i=0;i<=N;i++)
{
   xi = a + (dx*i);
   fn[i] = shift(x,f,xi,Ni,N,a,b,s);
}// tablichni stoinosti


for(i=0;i<=N;i++)
{
    xi = a + (dx*i);
//    printf("%g %g \n",xi,fn[i]);
}

for(i=0;i<=N;i++)
{
    f[i]=fn[i];
}

for(i=0;i<=N;i++)
{
   xi = a + (dx*i);
   printf("%g %g \n",xi,f[i]);
}


return 0;
}
