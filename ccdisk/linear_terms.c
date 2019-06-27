#include "parameters.h"
#include "variables.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void set_initial_values_linear_term(LinearTerm *lt,  double DQx,double DQy,double DQz,double DVQ)
{
    for (int ix=0;ix<=2*Nx; ix++)
	    for (int iy=0;iy<=2*Ny; iy++)
	        for (int iz=0;iz<=2*Nz; iz++)
                {
                double Qx=-Qx0+ix*DQx;
                double Qy=-Qy0+iy*DQy;
                double Q2 = Qx * Qx + Qy * Qy;
                double norm = 1 / sqrt(Q2);

                if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny)) == 0) 
                    continue;
                

                lt->x[ix][iy][iz] = vphi * -Qy * norm;
                lt->y[ix][iy][iz] = vphi * Qx * norm;
                lt->z[ix][iy][iz] = 0;
                }

    lt->x[Nx + 1][Ny][Nz] = 0.0;
    lt->y[Nx + 1][Ny][Nz] = vc;

    lt->x[Nx][Ny + 1][Nz] = -vc;
    lt->y[Nx][Ny + 1][Nz] = 0.0;

    lt->x[Nx - 1][Ny][Nz] = 0.0;
    lt->y[Nx - 1][Ny][Nz] = -vc;

    lt->x[Nx][Ny - 1][Nz] = vc;
    lt->y[Nx][Ny - 1][Nz] = 0.0;
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
						//  printf("Zero");
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

void to_be_translated()
{
            for (int ix = 0; ix <= Lx; ix++) {
            for (int iy = 0; iy <= Ly; iy++) {
                for (int iz = 0; iz <= Lz; iz++) {

                    double Qx = -Qx0 + ix * DQx;
                    double Qy = -Qy0 + iy * DQy;
                    double Qz = -Qz0 + iz * DQz;

                    double Q2 = Qx * Qx + Qy * Qy + Qz * Qz;
                    if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny) + (iz - Nz) * (iz - Nz)) == 0) {
                        continue;
                    }

                    if (abs(ix - Nx) == 1 && iy == Ny && iz == Nz) {
                        continue;
                    }

                    if (abs(iy - Ny) == 1 && ix == Nx && iz == Nz) {
                        continue;
                    }

                    double Q = sqrt(Q2);
                    double nx = Qx / Q;
                    double ny = Qy / Q;
                    double nz = Qz / Q;
                    double Qalpha = Qy * sin(alpha) + Qz * cos(alpha);

                    double nbx =  (bvxx[ix][iy][iz] - vbxx[ix][iy][iz]) * Qx + (bvxy[ix][iy][iz] - vbxy[ix][iy][iz]) * Qy +  (bvxz[ix][iy][iz] - vbxz[ix][iy][iz]) * Qz;
                    double nby =  (bvyx[ix][iy][iz] - vbyx[ix][iy][iz]) * Qx + (bvyy[ix][iy][iz] - vbyy[ix][iy][iz]) * Qy +  (bvyz[ix][iy][iz] - vbyz[ix][iy][iz]) * Qz;
                    double nbz =  (bvzx[ix][iy][iz] - vbzx[ix][iy][iz]) * Qx + (bvzy[ix][iy][iz] - vbzy[ix][iy][iz]) * Qy +  (bvzz[ix][iy][iz] - vbzz[ix][iy][iz]) * Qz;

                    double nvx =  (vvxx[ix][iy][iz] + bbxx[ix][iy][iz]) * Qx + (vvxy[ix][iy][iz] + bbxy[ix][iy][iz]) * Qy + (vvxz[ix][iy][iz] + bbxz[ix][iy][iz]) * Qz;
                    double nvy =  (vvyx[ix][iy][iz] + bbyx[ix][iy][iz]) * Qx + (vvyy[ix][iy][iz] + bbyy[ix][iy][iz]) * Qy + (vvyz[ix][iy][iz] + bbyz[ix][iy][iz]) * Qz;
                    double nvz =  (vvzx[ix][iy][iz] + bbzx[ix][iy][iz]) * Qx + (vvzy[ix][iy][iz] + bbzy[ix][iy][iz]) * Qy + (vvzz[ix][iy][iz] + bbzz[ix][iy][iz]) * Qz;

                                       
                    if ((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny) != 0)
                    {
                    double rho = sqrt(Qx*Qx+Qy*Qy);
                    double r = sqrt(Qx*Qx+ Qy*Qy+ Qz*Qz);

                    e.phy.x = -Qy/rho;
                    e.phy.y = Qx/rho;
                    e.phy.z = 0.0;

                    e.r.x = Qx/r;
                    e.r.y = Qy/r;
                    e.r.z = Qz/r;

                    e.theta.x = Qx*Qz/(rho*r);
                    e.theta.y = Qy*Qz/(rho*r);
                    e.theta.z = -(Qx*Qx+Qy*Qy)/(rho*r);
                    }
                    else
                    {
                    e.phy.x = 0.;
                    e.phy.y = 1.;
                    e.phy.z = 0.;

                    e.r.x = 0.;
                    e.r.y = 0.;
                    e.r.z = 1.;

                    e.theta.x = 1.;
                    e.theta.y = 0.;
                    e.theta.z = 0.;
                    }

                    double nbtheta = e.theta.x*nbx + e.theta.y*nby + e.theta.z*nbz;
                    double nbphy = e.phy.x*nbx + e.phy.y*nby + e.phy.z*nbz;
                    double nvtheta = e.theta.x*nvx + e.theta.y*nvy + e.theta.z*nvz;
                    double nvphy = e.phy.x*nvx + e.phy.y*nvy + e.phy.z*nvz;

                    double n4[3];
                    n4[0]=nbtheta;
                    n4[1]=nbphy;
                    n4[2]=nvtheta;
                    n4[3]=nvphy;
                                        
                    gsl_matrix * m = gsl_matrix_alloc (6, 6);
                    gsl_matrix * p = gsl_matrix_alloc (4, 6);
                    gsl_matrix_set_zero(m);
                    gsl_matrix_set_zero(p);

                    gsl_matrix_set (p, 0, 0, e.theta.x);
                    gsl_matrix_set (p, 0, 1, e.theta.y);
                    gsl_matrix_set (p, 0, 2, e.theta.z);
                    gsl_matrix_set (p, 0, 3, 0.0);
                    gsl_matrix_set (p, 0, 4, 0.0);
                    gsl_matrix_set (p, 0, 5, 0.0);

                    gsl_matrix_set (p, 1, 0, e.phy.x);
                    gsl_matrix_set (p, 1, 1, e.phy.y);
                    gsl_matrix_set (p, 1, 2, e.phy.z);
                    gsl_matrix_set (p, 1, 3, 0.0);
                    gsl_matrix_set (p, 1, 4, 0.0);
                    gsl_matrix_set (p, 1, 5, 0.0);

                    gsl_matrix_set (p, 2, 0, 0.0);
                    gsl_matrix_set (p, 2, 1, 0.0);
                    gsl_matrix_set (p, 2, 2, 0.0);
                    gsl_matrix_set (p, 2, 3, e.theta.x);
                    gsl_matrix_set (p, 2, 4, e.theta.y);
                    gsl_matrix_set (p, 2, 5, e.theta.z);

                    gsl_matrix_set (p, 3, 0, 0.0);
                    gsl_matrix_set (p, 3, 1, 0.0);
                    gsl_matrix_set (p, 3, 2, 0.0);
                    gsl_matrix_set (p, 3, 3, e.phy.x);
                    gsl_matrix_set (p, 3, 4, e.phy.y);
                    gsl_matrix_set (p, 3, 5, e.phy.z);

                    gsl_matrix_set (m, 0, 0, -nu_m * Q2);
                    gsl_matrix_set (m, 0, 1, 0);
                    gsl_matrix_set (m, 0, 2, 0);
                    gsl_matrix_set (m, 0, 3, -Qalpha);
                    gsl_matrix_set (m, 0, 4, 0);
                    gsl_matrix_set (m, 0, 5, 0);

                    gsl_matrix_set (m, 1, 0, 1);
                    gsl_matrix_set (m, 1, 1, -nu_m * Q2);
                    gsl_matrix_set (m, 1, 2, 0);
                    gsl_matrix_set (m, 1, 3, 0);
                    gsl_matrix_set (m, 1, 4, -Qalpha);
                    gsl_matrix_set (m, 1, 5, 0);


                    gsl_matrix_set (m, 2, 0, 0);
                    gsl_matrix_set (m, 2, 1, 0);
                    gsl_matrix_set (m, 2, 2, -nu_m * Q2);
                    gsl_matrix_set (m, 2, 3, 0);
                    gsl_matrix_set (m, 2, 4, 0);
                    gsl_matrix_set (m, 2, 5, -Qalpha);


                    gsl_matrix_set (m, 3, 0, Qalpha);
                    gsl_matrix_set (m, 3, 1, 0);
                    gsl_matrix_set (m, 3, 2, 0);
                    gsl_matrix_set (m, 3, 3, 2*ny*nx*(w+1) - nu_k*Q2);
                    gsl_matrix_set (m, 3, 4, -2*nx*nx*w + 2*w);
                    gsl_matrix_set (m, 3, 5, 0);


                    gsl_matrix_set (m, 4, 0, 0);
                    gsl_matrix_set (m, 4, 1, Qalpha);
                    gsl_matrix_set (m, 4, 2, 0);
                    gsl_matrix_set (m, 4, 3, 2*ny*ny*(w+1) - 2*w-1);
                    gsl_matrix_set (m, 4, 4, -2*nx*ny*w -nu_k*Q2);
                    gsl_matrix_set (m, 4, 5, 0);


                    gsl_matrix_set (m, 5, 0, 0);
                    gsl_matrix_set (m, 5, 1, 0);
                    gsl_matrix_set (m, 5, 2, Qalpha);
                    gsl_matrix_set (m, 5, 3, 2*ny*nz*(w+1));
                    gsl_matrix_set (m, 5, 4, -2*nx*nz*w);
                    gsl_matrix_set (m, 5, 5, -nu_k*Q2);

                    gsl_matrix * result = inverseOfMatrixMultiply(m,p);

                    double psi[3];
                    for (int ii = 0; ii <= 3; ii++)
                    {
                        for (int jj = 0; jj <= 3; jj++)
                        {
                            psi[ii]=gsl_matrix_get (result, ii, jj)*n4[jj];
                        }                        
                    }
                    bx[ix][iy][iz] = psi[0];
                    by[ix][iy][iz] = psi[1];                   
                    vx[ix][iy][iz] = psi[2];
                    vy[ix][iy][iz] = psi[3];
                    
                    }
            }
        }

}