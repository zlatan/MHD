#include <stdio.h>
#include <math.h>

int main(void)
{

double w=-2.0/3.0;
double Qx,Qz,Q2,Qa2;
double theta=0.0;
double b,c;
double x;

 for(Qx=-3.14;Qx<=3.14;Qx=Qx+0.1)
  for(Qz=-2.14;Qz<=2.14;Qz=Qz+0.1)
  {
 Qa2 = Qz*Qz*cos(theta)*cos(theta);
 Q2  = Qx*Qx + Qz*Qz;
   b = Qa2 + (1 +2*w)*( (Qx*Qx)/Q2 + 1 )*w;
   c = 2*Qa2*( (Qx*Qx)/Q2 + 1 )*w + Qa2*Qa2;
   x = -b + sqrt(b*b - c);
   printf("%g %g %g\n",Qx,Qz,x);
  }



return 0;
}
