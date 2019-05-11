#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

const double vc = 0.1;
const double vphi = 0.2;
const double alpha = 0.0;
const double w = -2. / 3.;
const double nu_k = 0.001;
const double nu_m = 0.001;


gsl_matrix * matrixMultiply(gsl_matrix *m, gsl_matrix *p)
{
gsl_matrix * tmp = gsl_matrix_alloc (4, 6);
gsl_matrix * result = gsl_matrix_alloc (4, 4);
gsl_matrix_set_zero(tmp);
gsl_matrix_set_zero(result);
gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, p, m, 0.0, tmp);
gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, tmp, p, 0.0, result);
gsl_matrix_free (tmp);
return result;
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

typedef struct vector
{
      double x;
      double y;
      double z;
} vector;

typedef struct unit_vector
{
      vector phy;
      vector r;
      vector theta;
} unit_vector;


int main() {
    double Qx0 = 3.0;
    double Qy0 = 3.0;
    double Qz0 = 3.0;

    int Nx = 6;
    int Ny = 6;
    int Nz = 6;

    int Nt = 100;
    double dt = 0.0025;
    
    int NDimx = 2 * Nx + 4, NDimy = 2 * Ny + 4, NDimz = 2 * Nz + 4;
    
    double vx[NDimx][NDimy][NDimz];
    double vy[NDimx][NDimy][NDimz];
    double vz[NDimx][NDimy][NDimz];

    double bx[NDimx][NDimy][NDimz];
    double by[NDimx][NDimy][NDimz];
    double bz[NDimx][NDimy][NDimz];
    
    unit_vector e;

    int Lx = 2 * Nx;
    int Ly = 2 * Ny;
    int Lz = 2 * Nz;

    double DQx = Qx0 / Nx;
    double DQy = Qy0 / Ny;
    double DQz = Qz0 / Nz;

    double DVQ = DQx * DQy * DQz / ((2 * M_PI) * (2 * M_PI) * (2 * M_PI));

    double bbyx[NDimx][NDimy][NDimz], bbyy[NDimx][NDimy][NDimz], bbyz[NDimx][NDimy][NDimz],
            bbxx[NDimx][NDimy][NDimz], bbxy[NDimx][NDimy][NDimz], bbxz[NDimx][NDimy][NDimz],
            bbzx[NDimx][NDimy][NDimz], bbzy[NDimx][NDimy][NDimz], bbzz[NDimx][NDimy][NDimz],
            vvxx[NDimx][NDimy][NDimz], vvxy[NDimx][NDimy][NDimz], vvxz[NDimx][NDimy][NDimz],
            vvyx[NDimx][NDimy][NDimz], vvyy[NDimx][NDimy][NDimz], vvyz[NDimx][NDimy][NDimz],
            vvzx[NDimx][NDimy][NDimz], vvzy[NDimx][NDimy][NDimz], vvzz[NDimx][NDimy][NDimz],
            bvxx[NDimx][NDimy][NDimz], bvxy[NDimx][NDimy][NDimz], bvxz[NDimx][NDimy][NDimz],
            bvyx[NDimx][NDimy][NDimz], bvyy[NDimx][NDimy][NDimz], bvyz[NDimx][NDimy][NDimz],
            bvzx[NDimx][NDimy][NDimz], bvzy[NDimx][NDimy][NDimz], bvzz[NDimx][NDimy][NDimz],
            vbxx[NDimx][NDimy][NDimz], vbxy[NDimx][NDimy][NDimz], vbxz[NDimx][NDimy][NDimz],
            vbyx[NDimx][NDimy][NDimz], vbyy[NDimx][NDimy][NDimz], vbyz[NDimx][NDimy][NDimz],
            vbzx[NDimx][NDimy][NDimz], vbzy[NDimx][NDimy][NDimz], vbzz[NDimx][NDimy][NDimz];

    for (int ix = 0; ix <= Lx; ix++) {
        for (int iy = 0; iy <= Ly; iy++) {
            for (int iz = 0; iz <= Nz; iz++) {
                double Qx = -Qx0 + ix * DQx;
                double Qy = -Qy0 + iy * DQy;

                double Q2 = Qx * Qx + Qy * Qy;
                double norm = 1 / sqrt(Q2);


                if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny)) == 0) {
                    continue;
                }

                vx[ix][iy][iz] = vphi * -Qy * norm;
                vy[ix][iy][iz] = vphi * Qx * norm;
                vz[ix][iy][iz] = 0;
            }
        }
    }

    vx[Nx + 1][Ny][Nz] = 0.0;
    vy[Nx + 1][Ny][Nz] = vc;


    vx[Nx][Ny + 1][Nz] = -vc;
    vy[Nx][Ny + 1][Nz] = 0.0;


    vx[Nx - 1][Ny][Nz] = 0.0;
    vy[Nx - 1][Ny][Nz] = -vc;

    vx[Nx][Ny - 1][Nz] = vc;
    vy[Nx][Ny - 1][Nz] = 0.0;


    for (int it = 1; it <= Nt; it++) {

        for (int ix = 0; ix <= Lx; ix++) {
            for (int iy = 0; iy <= Ly; iy++) {
                for (int iz = 0; iz <= Lz; iz++) {
                    bbxx[ix][iy][iz] = 0.0;
                    bbxy[ix][iy][iz] = 0.0;
                    bbxz[ix][iy][iz] = 0.0;
                    bbyx[ix][iy][iz] = 0.0;
                    bbyy[ix][iy][iz] = 0.0;
                    bbyz[ix][iy][iz] = 0.0;
                    bbzx[ix][iy][iz] = 0.0;
                    bbzy[ix][iy][iz] = 0.0;
                    bbzz[ix][iy][iz] = 0.0;

                    vvxx[ix][iy][iz] = 0.0;
                    vvxy[ix][iy][iz] = 0.0;
                    vvxz[ix][iy][iz] = 0.0;
                    vvyx[ix][iy][iz] = 0.0;
                    vvyy[ix][iy][iz] = 0.0;
                    vvyz[ix][iy][iz] = 0.0;
                    vvzx[ix][iy][iz] = 0.0;
                    vvzy[ix][iy][iz] = 0.0;
                    vvzz[ix][iy][iz] = 0.0;

                    bvxx[ix][iy][iz] = 0.0;
                    bvxy[ix][iy][iz] = 0.0;
                    bvxz[ix][iy][iz] = 0.0;
                    bvyx[ix][iy][iz] = 0.0;
                    bvyy[ix][iy][iz] = 0.0;
                    bvyz[ix][iy][iz] = 0.0;
                    bvzx[ix][iy][iz] = 0.0;
                    bvzy[ix][iy][iz] = 0.0;
                    bvzz[ix][iy][iz] = 0.0;

                    vbxx[ix][iy][iz] = 0.0;
                    vbxy[ix][iy][iz] = 0.0;
                    vbxz[ix][iy][iz] = 0.0;
                    vbyx[ix][iy][iz] = 0.0;
                    vbyy[ix][iy][iz] = 0.0;
                    vbyz[ix][iy][iz] = 0.0;
                    vbzx[ix][iy][iz] = 0.0;
                    vbzy[ix][iy][iz] = 0.0;
                    vbzz[ix][iy][iz] = 0.0;


                    #pragma omp for
                    for (int jx = 0; jx <= Lx; jx++) {
                        int kx = ix - jx + Nx;
                        if (kx < 0) continue;
                        if (kx > Lx) continue;
                        for (int jy = 0; jy <= Ly; jy++) {
                            int ky = iy - jy + Ny;
                            if (ky < 0) continue;
                            if (ky > Ly) continue;
                            for (int jz = 0; jz <= Lz; jz++) {
                                int kz = iz - jz + Nz;
                                if (kz < 0) continue;
                                if (kz > Lz) continue;

                                if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny) + (iz - Nz) * (iz - Nz)) == 0) {
                                    continue;
                                }
                                //bb
                                bbxx[ix][iy][iz] += bx[jx][jy][jz] * bx[kx][ky][kz] * DVQ;
                                bbxy[ix][iy][iz] += bx[jx][jy][jz] * by[kx][ky][kz] * DVQ;
                                bbxz[ix][iy][iz] += bx[jx][jy][jz] * bz[kx][ky][kz] * DVQ;
                                bbyx[ix][iy][iz] += by[jx][jy][jz] * bx[kx][ky][kz] * DVQ;
                                bbyy[ix][iy][iz] += by[jx][jy][jz] * by[kx][ky][kz] * DVQ;
                                bbyz[ix][iy][iz] += by[jx][jy][jz] * bz[kx][ky][kz] * DVQ;
                                bbzx[ix][iy][iz] += bz[jx][jy][jz] * bx[kx][ky][kz] * DVQ;
                                bbzy[ix][iy][iz] += bz[jx][jy][jz] * by[kx][ky][kz] * DVQ;
                                bbzz[ix][iy][iz] += bz[jx][jy][jz] * bz[kx][ky][kz] * DVQ;
                                //vv
                                vvxx[ix][iy][iz] += vx[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                vvxy[ix][iy][iz] += vx[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                vvxz[ix][iy][iz] += vx[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
                                vvyx[ix][iy][iz] += vy[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                vvyy[ix][iy][iz] += vy[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                vvyz[ix][iy][iz] += vy[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
                                vvzx[ix][iy][iz] += vz[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                vvzy[ix][iy][iz] += vz[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                vvzz[ix][iy][iz] += vz[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
                                //bv
                                bvxx[ix][iy][iz] += bx[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                bvxy[ix][iy][iz] += bx[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                bvxz[ix][iy][iz] += bx[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
                                bvyx[ix][iy][iz] += by[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                bvyy[ix][iy][iz] += by[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                bvyz[ix][iy][iz] += by[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
                                bvzx[ix][iy][iz] += bz[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                bvzy[ix][iy][iz] += bz[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                bvzz[ix][iy][iz] += bz[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
                                //vb
                                vbxx[ix][iy][iz] += vx[jx][jy][jz] * bx[kx][ky][kz] * DVQ;
                                vbxy[ix][iy][iz] += vx[jx][jy][jz] * by[kx][ky][kz] * DVQ;
                                vbxz[ix][iy][iz] += vx[jx][jy][jz] * bz[kx][ky][kz] * DVQ;
                                vbyx[ix][iy][iz] += vy[jx][jy][jz] * bx[kx][ky][kz] * DVQ;
                                vbyy[ix][iy][iz] += vy[jx][jy][jz] * by[kx][ky][kz] * DVQ;
                                vbyz[ix][iy][iz] += vy[jx][jy][jz] * bz[kx][ky][kz] * DVQ;
                                vbzx[ix][iy][iz] += vz[jx][jy][jz] * bx[kx][ky][kz] * DVQ;
                                vbzy[ix][iy][iz] += vz[jx][jy][jz] * by[kx][ky][kz] * DVQ;
                                vbzz[ix][iy][iz] += vz[jx][jy][jz] * bz[kx][ky][kz] * DVQ;
                            }
                        }
                    }
                }
            }
        }


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

        double EE = 0.0;
        //double vphi = 0.0;
        for (int ix = 0; ix <= Lx; ix++) {
            for (int iy = 0; iy <= Ly; iy++) {
                for (int iz = 0; iz <= Lz; iz++) {
                    double Qx = -Qx0 + ix * DQx;
                    double Qy = -Qy0 + iy * DQy;
                    EE += (
                        vx[ix][iy][iz] * vx[ix][iy][iz] + vy[ix][iy][iz] * vy[ix][iy][iz] + vz[ix][iy][iz] * vz[ix][iy][iz] 
                        +
                        bx[ix][iy][iz] * bx[ix][iy][iz] + by[ix][iy][iz] * by[ix][iy][iz] + bz[ix][iy][iz] * bz[ix][iy][iz]
                         )*DVQ;

                }
            }
        }
            printf("%i %.19f\n", it, EE);
    }
    return 0;


}