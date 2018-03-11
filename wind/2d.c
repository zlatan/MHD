#include <stdio.h>
#include <math.h>

int main(void)
{

int iy,ix;
int jy,jx;
int ky,kx;
double Ly=3.,Lx=3.;
double DQy,DQx,Q2,Q;
double Qy,Qx,Ky,Kx,K2,K;
int my=6,mx=6;
int Ny,Ny1,Nx,Nx1;
int i,k, N=3;
double f[5], g[5];
double x[5],xx;

Ny=2*my+1; Ny1=Ny+5;
Nx=2*mx+1; Nx1=Nx+5;
DQy=2*Ly/(Ny-1);
DQx=2*Lx/(Nx-1);

double Dt=.01,t;
double bx[Nx1][Ny1], by[Nx1][Ny1], ex[Nx1][Ny1], ey[Nx1][Ny1], vx[Nx1][Ny1], vy[Nx1][Ny1], ax[Nx1][Ny1], ay[Nx1][Ny1];
double qx[Nx1], Vx[Nx1][Ny1], Vy[Nx1][Ny1], Bx[Nx1][Ny1], By[Nx1][Ny1];
double nx, ny, nxp, nyp;
double lvx[Nx1][Ny1], lvy[Nx1][Ny1], lbx[Nx1][Ny1], lby[Nx1][Ny1], vn, Qb, Qpb, Qv,bb, nlv;// da se premesti
double energy;
double wc=0;//-2./3.;
double nuk=1/1000., num=1/137., Pk, Pm;
double v2,b2;

for(iy=1;iy<=Ny+4;iy++)
  for(ix=1;ix<=Nx+4;ix++)
 {
  vx[ix][iy]=4.3*10E-5;
  vy[ix][iy]=4.3*10E-5;
  bx[ix][iy]=4.3*10E-5;
  by[ix][iy]=4.3*10E-5;
  }



for(t=0;t<=5.5;t=t+Dt)
{

energy=0.;
Pk=nuk;
Pm=0.;

for(ix=1;ix<=Nx;ix++)
 {
  Qx=-Lx+(ix-1)*DQx; qx[ix]=Qx;
for(iy=1;iy<=Ny;iy++)
{
  Qy=-Ly+(iy-1)*DQy;
  //printf("Qy=%f iy=%i\n",Qy,iy);
// linear terms
  Q2=Qx*Qx+Qy*Qy; Q=sqrt(Q2); 
  if(Q != 0.) {nx=Qx/Q; ny=Qy/Q;} else {nx=0; ny=0;}
  ex[ix][iy]=ey[ix][iy]=ax[ix][iy]=ay[ix][iy]=0;
  ex[ix][iy]+= -num*Q2*bx[ix][iy]; 
  ey[ix][iy]+= +bx[ix][iy]-num*Q2*by[ix][iy];
  ax[ix][iy]+=              2*nx*ny*vx[ix][iy] +2*wc*nx*(ny*vx[ix][iy]-nx*vy[ix][iy]) +2*wc*vy[ix][iy] -nuk*Q2*vx[ix][iy];
  ay[ix][iy]+= -vx[ix][iy] +2*ny*ny*vx[ix][iy] +2*wc*ny*(ny*vx[ix][iy]-nx*vy[ix][iy]) -2*wc*vx[ix][iy]-nuk*Q2*vy[ix][iy];
 
//printf("ix=%i iy=%i ax=%g\n",ix,iy,ax[ix][iy]);
 lvx[ix][iy]=lvy[ix][iy]=lbx[ix][iy]=lby[ix][iy]=0;
for(jx=1;jx<=Nx;jx++)
   {
   Kx=-Lx+(jx-1)*DQx;
   kx=ix-jx+mx+1;
   if(kx<0) goto jump;
   if(kx>Nx)  goto jump;
 //  if(kx<0 && kx>Nx) continue; 
 //  if (Kx < -Lx && Kx > Lx ) continue; 
 for(jy=1;jy<=Ny;jy++)
   {
   Ky=-Ly+(jy-1)*DQy;
   ky=iy-jy+my+1;
   if(ky<0) goto jump; 
   if(ky>Ny) goto jump; 
  // if(ky<0 && ky>Ny) continue; 
 //  if (Ky < -Ly && Ky > Ly ) continue; 
// nonlinear terms
   
   vn=   vx[jx][jy]*Qx + vy[jx][jy]*Qy;
   Qb=   Qx*bx[kx][ky] + Qy*by[kx][ky];
   Qv=   Qx*vx[jx][jy] + Qy*vy[jx][jy];
   Qpb=  Kx*bx[kx][ky] + Ky*by[kx][ky];
   bb=   bx[jx][jy]*bx[kx][ky] +by[jx][jy]*by[kx][ky];
   lvx[ix][iy]+= vx[kx][ky]*vn + bx[jx][jy]*Qpb - Kx*bb;
   lvy[ix][iy]+= vy[kx][ky]*vn + by[jx][jy]*Qpb - Ky*bb;
   lbx[ix][iy]+= -vx[kx][ky]*Qb+ bx[kx][ky]*Qv;
   lby[ix][iy]+= -vy[kx][ky]*Qb+ by[kx][ky]*Qv;

jump: ;
 //printf("ix=%i iy=%i lvx=%g\n",ix,iy,lvx[ix][iy]); 
//printf("ix=%i iy=%i vn=%g\n",ix,iy,Qx);
    }// next jx
   }// next jy

if(Q != 0.) 
 {nlv=nx*lvx[Nx1][Ny1]+ny*lvy[Nx1][Ny1];
 lvx[Nx1][Ny1]=nx*nlv;
 lvy[Nx1][Ny1]=ny*nlv;}


vx[ix][iy]+=(ax[ix][iy]+lvx[ix][iy])*Dt; vy[ix][iy]+=(ay[ix][iy]+lvy[ix][iy])*Dt;
bx[ix][iy]+=(ex[ix][iy]+lbx[ix][iy])*Dt; by[ix][iy]+=(ey[ix][iy]+lby[ix][iy])*Dt;
v2=vx[ix][iy]*vx[ix][iy] + vy[ix][iy]*vy[ix][iy];
b2=by[ix][iy]*by[ix][iy] + bx[ix][iy]*bx[ix][iy];

xx=qx[ix]+Qy*Dt;

if(ix==1) 
 {
 x[0]=qx[1]; x[1]=qx[2]; x[2]=qx[3]; x[3]=qx[4];
 f[0]=vx[1][iy]; f[1]=vx[2][iy]; f[2]=vx[3][iy]; f[3]=vx[4][iy];
 }


if(ix>1 && ix<Nx-1) 
 {
 x[0]=qx[ix-1]; x[1]=qx[ix]; x[2]=qx[ix+1]; x[3]=qx[ix+2];
 f[0]=vx[ix-1][iy]; f[1]=vx[ix][iy]; f[2]=vx[ix+1][iy]; f[3]=vx[ix+2][iy];
 }

if(ix==Nx-1) 
 {
 x[0]=qx[Nx-1]; x[1]=qx[Nx]; x[2]=qx[Nx-1]; x[3]=qx[Nx-2];
 f[0]=vx[Nx-1][iy]; f[1]=vx[Nx][iy]; f[2]=vx[Nx-1][iy]; f[3]=vx[Nx-2][iy];
 }

if(ix==Nx) 
 {
 x[0]=qx[Nx]; x[1]=qx[Nx-1]; x[2]=qx[Nx-2]; x[3]=qx[Nx-3];
 f[0]=vx[Nx][iy]; f[1]=vx[Nx-1][iy]; f[2]=vx[Nx-2][iy]; f[3]=vx[Nx-3][iy];
 }

/*
for(i=1;i<=N;i++)
 printf("%g\n",f[i]);
*/

for(k=0;k<N; k++)
for(i=k+1;i<=N;i++) 
{
f[i]=( (x[k]-xx)*f[i] - (x[i]-xx)*f[k] )/( (x[k]-xx)- (x[i]-xx) );
}
Vx[ix][iy]=f[N]; // srestu vjatyra

///////////////////////////////////////////////////////////////////////////////////

if(ix==1) 
 {
 x[0]=qx[1]; x[1]=qx[2]; x[2]=qx[3]; x[3]=qx[4];
 f[0]=bx[1][iy]; f[1]=bx[2][iy]; f[2]=bx[3][iy]; f[3]=bx[4][iy];
 }


if(ix>1 && ix<Nx-1) 
 {
 x[0]=qx[ix-1]; x[1]=qx[ix]; x[2]=qx[ix+1]; x[3]=qx[ix+2];
 f[0]=bx[ix-1][iy]; f[1]=bx[ix][iy]; f[2]=bx[ix+1][iy]; f[3]=bx[ix+2][iy];
 }

if(ix==Nx-1) 
 {
 x[0]=qx[Nx-1]; x[1]=qx[Nx]; x[2]=qx[Nx-1]; x[3]=qx[Nx-2];
 f[0]=bx[Nx-1][iy]; f[1]=bx[Nx][iy]; f[2]=bx[Nx-1][iy]; f[3]=bx[Nx-2][iy];
 }

if(ix==Nx) 
 {
 x[0]=qx[Nx]; x[1]=qx[Nx-1]; x[2]=qx[Nx-2]; x[3]=qx[Nx-3];
 f[0]=bx[Nx][iy]; f[1]=bx[Nx-1][iy]; f[2]=bx[Nx-2][iy]; f[3]=bx[Nx-3][iy];
 }


for(k=0;k<N; k++)
for(i=k+1;i<=N;i++) 
{
f[i]=( (x[k]-xx)*f[i] - (x[i]-xx)*f[k] )/( (x[k]-xx)- (x[i]-xx) );
}
Bx[ix][iy]=f[N]; // srestu vjatyra

//////////////////////////////////////////////////////////////////////////////


if(ix==1) 
 {
 x[0]=qx[1]; x[1]=qx[2]; x[2]=qx[3]; x[3]=qx[4];
 f[0]=vy[1][iy]; f[1]=vy[2][iy]; f[2]=vy[3][iy]; f[3]=vy[4][iy];
 }


if(ix>1 && ix<Nx-1) 
 {
 x[0]=qx[ix-1]; x[1]=qx[ix]; x[2]=qx[ix+1]; x[3]=qx[ix+2];
 f[0]=vy[ix-1][iy]; f[1]=vy[ix][iy]; f[2]=vy[ix+1][iy]; f[3]=vy[ix+2][iy];
 }

if(ix==Nx-1) 
 {
 x[0]=qx[Nx-1]; x[1]=qx[Nx]; x[2]=qx[Nx-1]; x[3]=qx[Nx-2];
 f[0]=vy[Nx-1][iy]; f[1]=vy[Nx][iy]; f[2]=vy[Nx-1][iy]; f[3]=vy[Nx-2][iy];
 }

if(ix==Nx) 
 {
 x[0]=qx[Nx]; x[1]=qx[Nx-1]; x[2]=qx[Nx-2]; x[3]=qx[Nx-3];
 f[0]=vy[Nx][iy]; f[1]=vy[Nx-1][iy]; f[2]=vy[Nx-2][iy]; f[3]=vy[Nx-3][iy];
 }



for(k=0;k<N; k++)
for(i=k+1;i<=N;i++) 
{
f[i]=( (x[k]-xx)*f[i] - (x[i]-xx)*f[k] )/( (x[k]-xx)- (x[i]-xx) );
}
Vy[ix][iy]=f[N]; // srestu vjatyra

////////////////////////////////////////////////////////


if(ix==1) 
 {
 x[0]=qx[1]; x[1]=qx[2]; x[2]=qx[3]; x[3]=qx[4];
 f[0]=by[1][iy]; f[1]=by[2][iy]; f[2]=by[3][iy]; f[3]=by[4][iy];
 }


if(ix>1 && ix<Nx-1) 
 {
 x[0]=qx[ix-1]; x[1]=qx[ix]; x[2]=qx[ix+1]; x[3]=qx[ix+2];
 f[0]=by[ix-1][iy]; f[1]=by[ix][iy]; f[2]=by[ix+1][iy]; f[3]=by[ix+2][iy];
 }

if(ix==Nx-1) 
 {
 x[0]=qx[Nx-1]; x[1]=qx[Nx]; x[2]=qx[Nx-1]; x[3]=qx[Nx-2];
 f[0]=by[Nx-1][iy]; f[1]=by[Nx][iy]; f[2]=by[Nx-1][iy]; f[3]=by[Nx-2][iy];
 }

if(ix==Nx) 
 {
 x[0]=qx[Nx]; x[1]=qx[Nx-1]; x[2]=qx[Nx-2]; x[3]=qx[Nx-3];
 f[0]=by[Nx][iy]; f[1]=by[Nx-1][iy]; f[2]=by[Nx-2][iy]; f[3]=by[Nx-3][iy];
 }



for(k=0;k<N; k++)
for(i=k+1;i<=N;i++) 
{
f[i]=( (x[k]-xx)*f[i] - (x[i]-xx)*f[k] )/( (x[k]-xx)- (x[i]-xx) );
}
By[ix][iy]=f[N]; // srestu vjatyra

/////////////////////////////////////////////////////

energy+= (v2 + b2)/2.; 


 } // iy
} // ix


for(ix=1;ix<=Nx;ix++)
 {
 for(iy=1;iy<=Ny;iy++) 
  {
  vx[ix][iy]=Vx[ix][iy]; 
  vy[ix][iy]=Vy[ix][iy];
  bx[ix][iy]=Bx[ix][iy];
  by[ix][iy]=By[ix][iy];
  }// iy
 }// ix otchetohme vyatyra 

printf("%g %g\n",t,log(energy));
}
return 0;
}
