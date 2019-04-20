#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


const double nu = 0.003;
const double vc = 0.1;
const double vphi = 0.2;


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
    
    unit_vector e;
    Matrixp pmatrix;

    int Lx = 2 * Nx;
    int Ly = 2 * Ny;
    int Lz = 2 * Nz;

    double DQx = Qx0 / Nx;
    double DQy = Qy0 / Ny;
    double DQz = Qz0 / Nz;

    double DVQ = DQx * DQy * DQz / ((2 * M_PI) * (2 * M_PI) * (2 * M_PI));

    double vvxx[NDimx][NDimy][NDimz], vvxy[NDimx][NDimy][NDimz], vvxz[NDimx][NDimy][NDimz],
           vvyx[NDimx][NDimy][NDimz], vvyy[NDimx][NDimy][NDimz], vvyz[NDimx][NDimy][NDimz],
           vvzx[NDimx][NDimy][NDimz], vvzy[NDimx][NDimy][NDimz], vvzz[NDimx][NDimy][NDimz];

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
                    vvxx[ix][iy][iz] = 0.0;
                    vvxy[ix][iy][iz] = 0.0;
                    vvxz[ix][iy][iz] = 0.0;
                    vvyx[ix][iy][iz] = 0.0;
                    vvyy[ix][iy][iz] = 0.0;
                    vvyz[ix][iy][iz] = 0.0;
                    vvzx[ix][iy][iz] = 0.0;
                    vvzy[ix][iy][iz] = 0.0;
                    vvzz[ix][iy][iz] = 0.0;


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

                                vvxx[ix][iy][iz] += vx[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                vvxy[ix][iy][iz] += vx[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                vvxz[ix][iy][iz] += vx[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
                                vvyx[ix][iy][iz] += vy[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                vvyy[ix][iy][iz] += vy[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                vvyz[ix][iy][iz] += vy[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
                                vvzx[ix][iy][iz] += vz[jx][jy][jz] * vx[kx][ky][kz] * DVQ;
                                vvzy[ix][iy][iz] += vz[jx][jy][jz] * vy[kx][ky][kz] * DVQ;
                                vvzz[ix][iy][iz] += vz[jx][jy][jz] * vz[kx][ky][kz] * DVQ;
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

                    double fnvx = vvxx[ix][iy][iz] * Qx + vvxy[ix][iy][iz] * Qy + vvxz[ix][iy][iz] * Qz;
                    double fnvy = vvyx[ix][iy][iz] * Qx + vvyy[ix][iy][iz] * Qy + vvyz[ix][iy][iz] * Qz;
                    double fnvz = vvzx[ix][iy][iz] * Qx + vvzy[ix][iy][iz] * Qy + vvzz[ix][iy][iz] * Qz;

                    vx[ix][iy][iz] = fnvx/(nu*Q2);
                    vy[ix][iy][iz] = fnvy/(nu*Q2);
                    vz[ix][iy][iz] = fnvz/(nu*Q2);

                    
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

                    gsl_matrix * result = matrixMultiply(m,p);


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
                    EE += (vx[ix][iy][iz] * vx[ix][iy][iz] + vy[ix][iy][iz] * vy[ix][iy][iz] +
                          vz[ix][iy][iz] * vz[ix][iy][iz])*DVQ;

                }
            }
        }
            printf("%i %g\n", it, EE);
    }
    return 0;


}