#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>


const double nu = 0.0001;
const double vc = 2.5;
const double vphi = 0.2;

int main() {
    double Qx0 = 4.0;
    double Qy0 = 4.0;
    double Qz0 = 4.0;

    int Nx = 6;
    int Ny = 6;
    int Nz = 6;

    int Nt = 10000;
    double dt = 0.0025;
    
    int NDimx = 2 * Nx + 4, NDimy = 2 * Ny + 4, NDimz = 2 * Nz + 4;
    
    double vx[NDimx][NDimy][NDimz];
    double vy[NDimx][NDimy][NDimz];
    double vz[NDimx][NDimy][NDimz];
    
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

                    double flvx = -nu * Q2 * vx[ix][iy][iz];
                    double flvy = -nu * Q2 * vy[ix][iy][iz];
                    double flvz = -nu * Q2 * vz[ix][iy][iz];


                    double fvx = flvx + fnvx;
                    double fvy = flvy + fnvy;
                    double fvz = flvz + fnvz;

                    double fvn = fvx * nx + fvy * ny + fvz * nz;


                    fvx -= fvn * nx;
                    fvy -= fvn * ny;
                    fvz -= fvn * nz;


                    vx[ix][iy][iz] += fvx * dt;
                    vy[ix][iy][iz] += fvy * dt;
                    vz[ix][iy][iz] += fvz * dt;

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
//                    double Qz = -Qz0 + iz * DQz;
                    EE += vx[ix][iy][iz] * vx[ix][iy][iz] + vy[ix][iy][iz] * vy[ix][iy][iz] +
                          vz[ix][iy][iz] * vz[ix][iy][iz];

//                double vphi = (-Qy * vx[ix][iy][Nz] + Qx * vy[ix][iy][Nz]) / sqrt(Qx * Qx + Qy * Qy);
//                printf("%f %f %f\n", Qx, Qy, vphi);
                }
            }
        }
            printf("%i %g\n", it, EE);

    }

    return 0;
}
/*
 * vx[Nx + 1][Ny][Nz] = 0.0;
    vy[Nx + 1][Ny][Nz] = vc;


    vx[Nx][Ny + 1][Nz] = -vc;
    vy[Nx][Ny + 1][Nz] = 0.0;


    vx[Nx - 1][Ny][Nz] = 0.0;
    vy[Nx - 1][Ny][Nz] = -vc;

    vx[Nx][Ny - 1][Nz] = vc;
    vy[Nx][Ny - 1][Nz] = 0.0;

*/
