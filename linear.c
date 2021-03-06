#include "parameters.h"
#include "variables.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


void integrate_over_t(LinearTerm *b, LinearTerm *v, LinearTerm *fv, LinearTerm *fb, double dt)
{

	for(int ix=0;ix<=2*Nx; ix++)
		for(int iy=0;iy<=2*Ny; iy++)
			for(int iz=0;iz<=2*Nz; iz++)
			{
			v->x[ix][iy][iz]+=fv->x[ix][iy][iz]*dt;
 			v->y[ix][iy][iz]+=fv->y[ix][iy][iz]*dt;
			v->z[ix][iy][iz]+=fv->z[ix][iy][iz]*dt;
		        b->x[ix][iy][iz]+=fb->x[ix][iy][iz]*dt;
			b->y[ix][iy][iz]+=fb->y[ix][iy][iz]*dt;
			b->z[ix][iy][iz]+=fb->z[ix][iy][iz]*dt;
			}
}

void calculate_linear_term(NonlinearTerm *bb, NonlinearTerm *vv, NonlinearTerm *bv, NonlinearTerm *vb, LinearTerm *b, LinearTerm *v, LinearTerm *fv, LinearTerm *fb, double DQx,double DQy,double DQz,double DVQ)
{
	double betax,betaz;
	double fnvx,fnvy,fnvz,fnbx,fnby,fnbz;
	double Q2,Qbeta,Q;
	betax = sin(alpha);
	betaz = cos(alpha);


			for(int ix=0;ix<=2*Nx; ix++)
			{
			double Qx=-Qx0+ix*DQx;
				for(int iy=0;iy<=2*Nx; iy++)
				{
				double Qy=-Qy0+iy*DQy;
					for(int iz=0;iz<=2*Nx; iz++)
					{
					double Qz=-Qz0+iz*DQz;
   				        double rnx, rny, rnz;
					double fvn,fbn,aaa,fnvx,fnvy,fnvz,fnbx,fnby,fnbz,flvx,flvy,flvz,flbx,flby,flbz;

					 Q2=Qx*Qx+Qy*Qy+Qz*Qz;
					 if ( ((ix - Nx)*(ix - Nx) + (iy-Ny)*(iy-Ny) + (iz-Nz)*(iz-Nz)) == 0 )
					 {
						 printf("Zero");
						 continue;
					 }
					 Q=sqrt(Q2);
					//  printf("Q=%g it=>%i",Q,it);
					 rnx=Qx/Q; rny=Qy/Q; rnz=Qz/Q;
					 Qbeta=Qx*betax+Qz*betaz; // betay=0.

					 fnvx = (vv->xx[ix][iy][iz] + bb->xx[ix][iy][iz] )* Qx + (vv->xy[ix][iy][iz] + bb->xy[ix][iy][iz]) * Qy + (vv->xz[ix][iy][iz] + bb->xz[ix][iy][iz]) * Qz;
					 fnvy = (vv->yx[ix][iy][iz] + bb->yx[ix][iy][iz] )* Qx + (vv->yy[ix][iy][iz] + bb->yy[ix][iy][iz]) * Qy + (vv->yz[ix][iy][iz] + bb->yz[ix][iy][iz]) * Qz;
					 fnvz = (vv->zx[ix][iy][iz] + bb->zx[ix][iy][iz] )* Qx + (vv->zy[ix][iy][iz] + bb->zy[ix][iy][iz]) * Qy + (vv->zz[ix][iy][iz] + bb->zz[ix][iy][iz]) * Qz;

					//  printf("fnvx=%g fnvy=%g fnvz=%g iteration=>%i\n",fnvx,fnvy,fnvz,it);
					//  printf("rnx=%g rny=%g rnz=%g iteration=>%i\n",rnx,rny,rnz,it);

					 fnbx = (bv->xx[ix][iy][iz] - vb->xx[ix][iy][iz] )* Qx + (bv->xy[ix][iy][iz] - vb->xy[ix][iy][iz]) * Qy + (bv->xz[ix][iy][iz] - vb->xz[ix][iy][iz]) * Qz;
					 fnby = (bv->yx[ix][iy][iz] - vb->yx[ix][iy][iz] )* Qx + (bv->yy[ix][iy][iz] - vb->yy[ix][iy][iz]) * Qy + (bv->yz[ix][iy][iz] - vb->yz[ix][iy][iz]) * Qz;
					 fnbz = (bv->zx[ix][iy][iz] - vb->zx[ix][iy][iz] )* Qx + (bv->zy[ix][iy][iz] - vb->zy[ix][iy][iz]) * Qy + (bv->zz[ix][iy][iz] - vb->zz[ix][iy][iz]) * Qz;

					 fvn = fnvx*rnx + fnvy*rny + fnvz*rnz;
					 fbn = fnbx*rnx + fnby*rny + fnbz*rnz;

					//  printf("fvn=%g it=>%i",fvn,it);

					 fnvx-= fvn*rnx;
					 fnvy-= fvn*rny;
					 fnvz-= fvn*rnz;

					//  printf("fnvx=%g fnvy=%g fnvz=%g iteration=>%i\n",fnvx,fnvy,fnvz,it);

					 fnbx-= fbn*rnx;
					 fnby-= fbn*rny;
					 fnbz-= fbn*rnz;

					 aaa  = 2.0*(rny*v->x[ix][iy][iz] - omega*(rnx*v->y[ix][iy][iz]-rny*v->x[ix][iy][iz]));
					 flvx =                   aaa*rnx + 2.0*omega*v->y[ix][iy][iz] + Qbeta*b->x[ix][iy][iz] - rnuk*Q2*v->x[ix][iy][iz];
					 flvy =-v->x[ix][iy][iz] + aaa*rny - 2.0*omega*v->x[ix][iy][iz] + Qbeta*b->y[ix][iy][iz] - rnuk*Q2*v->y[ix][iy][iz];
					 flvz =                   aaa*rnz                              + Qbeta*b->z[ix][iy][iz] - rnuk*Q2*v->z[ix][iy][iz];

					 flbx =      		  - Qbeta*v->x[ix][iy][iz] - rnum*Q2*v->x[ix][iy][iz];
					 flby = b->x[ix][iy][iz]  - Qbeta*v->y[ix][iy][iz] - rnum*Q2*v->y[ix][iy][iz];
					 flbz =    		  - Qbeta*v->z[ix][iy][iz] - rnum*Q2*v->z[ix][iy][iz];


					 fvn = flvx*rnx + flvy*rny + flvz*rnz;
					 fbn = flbx*rnx + flby*rny + flbz*rnz;

					//  printf("fvn=%g fbn=%g iteration=>%i\n",fvn,fbn,it);

					 flvx-= fvn*rnx;
					 flvy-= fvn*rny;
					 flvz-= fvn*rnz;

					 flbx-= fbn*rnx;
					 flby-= fbn*rny;
					 flbz-= fbn*rnz;

					//  printf("flvx=%g flvy=%g flvz=%g iteration=>%i\n",flvx,flvy,flvz,it);

					 fv->x[ix][iy][iz] = flvx+fnvx;
					 fv->y[ix][iy][iz] = flvy+fnvy;
					 fv->z[ix][iy][iz] = flvz+fnvz;

					//  printf("x->%g y->%g z->%g \n",fvx[ix][iy][iz],fvy[ix][iy][iz],fvy[ix][iy][iz]);

					 fb->x[ix][iy][iz] = flbx+fnbx;
					 fb->y[ix][iy][iz] = flby+fnby;
					 fb->z[ix][iy][iz] = flbz+fnbz;

					}
				}
			}
}
