#include "parameters.h"
#include "variables.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

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

gsl_matrix * inverseOfMatrixMultiply(gsl_matrix *m, gsl_matrix *p)
{
gsl_matrix * tmp = gsl_matrix_alloc (4, 6);
gsl_matrix * result = gsl_matrix_alloc (4, 4);
gsl_matrix * inverse = gsl_matrix_alloc (4, 4);
gsl_permutation * per = gsl_permutation_alloc (4);

int s;

gsl_matrix_set_zero(tmp);
gsl_matrix_set_zero(result);
gsl_matrix_set_zero(inverse);


gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, p, m, 0.0, tmp);
gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, tmp, p, 0.0, result);
gsl_matrix_free (tmp);

gsl_linalg_LU_decomp (result, per, &s);    
gsl_linalg_LU_invert (result, per, inverse);
return inverse;
}


void calculate_linear_term(NonlinearTerm *bb, NonlinearTerm *vv, NonlinearTerm *bv, NonlinearTerm *vb, LinearTerm *b, LinearTerm *v, double DQx,double DQy,double DQz,double DVQ)
{
    unit_vector e;

			for(int ix=0;ix<=2*Nx; ix++)
			{
			double Qx=-Qx0+ix*DQx;
				for(int iy=0;iy<=2*Nx; iy++)
				{
				double Qy=-Qy0+iy*DQy;
					for(int iz=0;iz<=2*Nx; iz++)
					{
					double Qz=-Qz0+iz*DQz;

					double Q2=Qx*Qx+Qy*Qy+Qz*Qz;
                    double Q = sqrt(Q2);
                    double Qalpha = Qy * sin(aleph) + Qz * cos(aleph);

                    double nx = Qx/Q;
                    double ny = Qy/Q;
                    double nz = Qz/Q;
					double nvx,nvy,nvz,nbx,nby,nbz;

                    if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny) + (iz - Nz) * (iz - Nz)) == 0) {
                        continue;
                    }

                    if (abs(ix - Nx) == 1 && iy == Ny && iz == Nz) {
                        continue;
                    }

                    if (abs(iy - Ny) == 1 && ix == Nx && iz == Nz) {
                        continue;
                    }


					 nvx = (vv->xx[ix][iy][iz] + bb->xx[ix][iy][iz] )* Qx + (vv->xy[ix][iy][iz] + bb->xy[ix][iy][iz]) * Qy + (vv->xz[ix][iy][iz] + bb->xz[ix][iy][iz]) * Qz;
					 nvy = (vv->yx[ix][iy][iz] + bb->yx[ix][iy][iz] )* Qx + (vv->yy[ix][iy][iz] + bb->yy[ix][iy][iz]) * Qy + (vv->yz[ix][iy][iz] + bb->yz[ix][iy][iz]) * Qz;
					 nvz = (vv->zx[ix][iy][iz] + bb->zx[ix][iy][iz] )* Qx + (vv->zy[ix][iy][iz] + bb->zy[ix][iy][iz]) * Qy + (vv->zz[ix][iy][iz] + bb->zz[ix][iy][iz]) * Qz;

					 nbx = (bv->xx[ix][iy][iz] - vb->xx[ix][iy][iz] )* Qx + (bv->xy[ix][iy][iz] - vb->xy[ix][iy][iz]) * Qy + (bv->xz[ix][iy][iz] - vb->xz[ix][iy][iz]) * Qz;
					 nby = (bv->yx[ix][iy][iz] - vb->yx[ix][iy][iz] )* Qx + (bv->yy[ix][iy][iz] - vb->yy[ix][iy][iz]) * Qy + (bv->yz[ix][iy][iz] - vb->yz[ix][iy][iz]) * Qz;
					 nbz = (bv->zx[ix][iy][iz] - vb->zx[ix][iy][iz] )* Qx + (bv->zy[ix][iy][iz] - vb->zy[ix][iy][iz]) * Qy + (bv->zz[ix][iy][iz] - vb->zz[ix][iy][iz]) * Qz;


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

                    double n4[4];
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

                    gsl_matrix_set (m, 0, 0, -rnum * Q2);
                    gsl_matrix_set (m, 0, 1, 0);
                    gsl_matrix_set (m, 0, 2, 0);
                    gsl_matrix_set (m, 0, 3, -Qalpha);
                    gsl_matrix_set (m, 0, 4, 0);
                    gsl_matrix_set (m, 0, 5, 0);

                    gsl_matrix_set (m, 1, 0, 1);
                    gsl_matrix_set (m, 1, 1, -rnum * Q2);
                    gsl_matrix_set (m, 1, 2, 0);
                    gsl_matrix_set (m, 1, 3, 0);
                    gsl_matrix_set (m, 1, 4, -Qalpha);
                    gsl_matrix_set (m, 1, 5, 0);


                    gsl_matrix_set (m, 2, 0, 0);
                    gsl_matrix_set (m, 2, 1, 0);
                    gsl_matrix_set (m, 2, 2, -rnum * Q2);
                    gsl_matrix_set (m, 2, 3, 0);
                    gsl_matrix_set (m, 2, 4, 0);
                    gsl_matrix_set (m, 2, 5, -Qalpha);


                    gsl_matrix_set (m, 3, 0, Qalpha);
                    gsl_matrix_set (m, 3, 1, 0);
                    gsl_matrix_set (m, 3, 2, 0);
                    gsl_matrix_set (m, 3, 3, 2*ny*nx*(omega+1) - rnuk*Q2);
                    gsl_matrix_set (m, 3, 4, -2*nx*nx*omega + 2*omega);
                    gsl_matrix_set (m, 3, 5, 0);


                    gsl_matrix_set (m, 4, 0, 0);
                    gsl_matrix_set (m, 4, 1, Qalpha);
                    gsl_matrix_set (m, 4, 2, 0);
                    gsl_matrix_set (m, 4, 3, 2*ny*ny*(omega+1) - 2*omega-1);
                    gsl_matrix_set (m, 4, 4, -2*nx*ny*omega -rnuk*Q2);
                    gsl_matrix_set (m, 4, 5, 0);


                    gsl_matrix_set (m, 5, 0, 0);
                    gsl_matrix_set (m, 5, 1, 0);
                    gsl_matrix_set (m, 5, 2, Qalpha);
                    gsl_matrix_set (m, 5, 3, 2*ny*nz*(omega+1));
                    gsl_matrix_set (m, 5, 4, -2*nx*nz*omega);
                    gsl_matrix_set (m, 5, 5, -rnuk*Q2);

                    gsl_matrix * result = inverseOfMatrixMultiply(m,p);

                    double psi[3];
                    for (int ii = 0; ii <= 3; ii++)
                    {
                        for (int jj = 0; jj <= 3; jj++)
                        {
                            psi[ii]=gsl_matrix_get (result, ii, jj)*n4[jj];
                        }                        
                    }
                    b->x[ix][iy][iz] = psi[0];
                    b->y[ix][iy][iz] = psi[1];                   
                    v->x[ix][iy][iz] = psi[2];
                    v->y[ix][iy][iz] = psi[3];

                }					
			}
		}   
}