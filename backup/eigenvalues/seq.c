#include <stdio.h>
#include <math.h>

int main(void)
{

double b,c,q,x;
double w=-2.0/3.0;


for(q=-1.2;q<=1.2;q=q+0.01)
{
b = q*q + (1+2*w)*w;
c = (q*q + 2*w)*q*q;
x = -b + sqrt(b*b-c);
printf("%g %g\n",q,x);
}

/*
b = (1+2*w)*w;
c = 2*w;
x = -b + sqrt(b*b-c);
printf("0 -> %g\n",x);
*/


return 0;
}
