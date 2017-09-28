#include <stdio.h>
#include <math.h> 


int D=4;

const double w=-2./3.,theta=M_PI/3;
const double K_y=0.05;

void force(double f[D],double y[D], double t,double K_y, double K_z)
{
double K2,K_x,K_a;
K2=(K_y*K_y*(1+t*t)+K_z*K_z);
K_x=-t*K_y;
K_a=K_y*sin(theta)+K_z*cos(theta);

f[1]= -K_a*y[2];				//b_x
f[2]= (2*K_y*K_x*y[2])/K2 + K_a*y[1]+2*w*y[4] + ((2*w)/K2)*(y[2]*K_x*K_y -y[4]*K_x*K_x);		//v_x
f[3]= y[1] -K_a*y[4];				//b_y
f[4]= (2*K_y*K_y*y[2])/K2 -y[2]+ K_a*y[3] -2*w*y[2] + ((2*w)/K2)*(y[2]*K_y*K_y -y[4]*K_x*K_y);	//v_y
}

double calc(double K_y, double K_z,double phi)
{

double a,b,y0[D],h2,k1[D],k2[D],k3[D],k4[D],dy[D],h,y[D],yc[D],t;
int N,n,d;
double tau;

a=-13/K_y;
b=-a;
N=10000;
h=(b-a)/(N-1);
h2=h/2;

y0[1]=(1/a)*(cos(K_y*a-phi));//b_x
y0[2]=-1*(1/a)*(sin(K_y*a-phi));//v_x
y0[3]=cos(K_y*a-phi);//b_y
y0[4]=-1*sin(K_y*a-phi);//v_y


t=a;

for(d=1;d<=D;d++) {y[d]=y0[d];}

for(n=1;n<N;n++)
{
t=a+h*(n-1);// t_n
force(k1,y,t,K_y,K_z);
for(d=1;d<=D;d++) yc[d]=y[d]+h2*k1[d];
force(k2,yc,t+h2,K_y,K_z);
for(d=1;d<=D;d++) yc[d]=y[d]+h2*k2[d];
force(k3,yc,t+h2,K_y,K_z);
for(d=1;d<=D;d++) yc[d]=y[d]+h*k3[d];
force(k4,yc,t+h,K_y,K_z);

for(d=1;d<=D;d++)  y[d]+=h*(k1[d]+k4[d]+(k2[d]+k3[d])*2)/6; 

tau=t+h;
}
//printf("phi=%f %f\n",phi,(y[1]*y[1]+y[2]*y[2]+y[3]*y[3]+y[4]*y[4]) );
return ( y[1]*y[1] + y[2]*y[2] + y[3]*y[3] + y[4]*y[4] );
}


int main(void)
{
double K_z,phi,Amp;
int i,Np=361;

  for(K_z=0;K_z<5;K_z+=0.1)
    {		
	for(i=1;i<Np;i++)
 	{
		phi=(M_PI*(i-1))/Np;
		Amp+=calc(K_y,K_z,phi);
	 }
    Amp=Amp/(Np-1);
 printf("%f %f\n",K_z,Amp);
   Amp=0;
   phi=0;
   }


return 0;
}



