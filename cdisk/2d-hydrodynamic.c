#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>


const double nu = 0.000012;


int main() {
    double Qx0 = 4.0;
    double Qy0 = 4.0;

    int Nx = 16;
    int Ny = 16;

    int Nt = 5000;
    double dt = 0.0025;

    int Ni = 4;

    int NDimx = 2 * Nx + 4, NDimy = 2 * Ny + 4;

    int NN = 2 * Nx;

    double vx[NDimx][NDimy];
    double vy[NDimx][NDimy];

    int Lx = 2 * Nx;
    int Ly = 2 * Ny;

    double DQx = Qx0 / Nx;
    double DQy = Qy0 / Ny;

    double DVQ = DQx * DQy / ((2 * M_PI) * (2 * M_PI));

    double vvxx[NDimx][NDimy], vvxy[NDimx][NDimy];
    double vvyx[NDimx][NDimy], vvyy[NDimx][NDimy];
    double vvzx[NDimx][NDimy], vvzy[NDimx][NDimy];


    for (int ix = 0; ix <= Lx; ix++) {
        for (int iy = 0; iy <= Ly; iy++) {
            double Qx = -Qx0 + ix * DQx;
            double Qy = -Qy0 + iy * DQy;

            vx[ix][iy] = arc4random_uniform(10) * 1E-2;
            vy[ix][iy] = arc4random_uniform(10) * 1E-2;


            double Q2 = Qx * Qx + Qy * Qy;
            double Q = sqrt(Q2);

            if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny)) == 0) {
                continue;
            }
            double nx = Qx / Q;
            double ny = Qy / Q;

            double vn = vx[ix][iy] * nx + vy[ix][iy] * ny;

            vx[ix][iy] -= vn * nx;
            vy[ix][iy] -= vn * ny;
        }
    }


    for (int it = 1; it <= Nt; it++) {
        vx[Nx+1][0] = 1;
        vy[Nx-1][0] = -1;

        for (int ix = 0; ix <= Lx; ix++) {
            for (int iy = 0; iy <= Ly; iy++) {
                vvxx[ix][iy] = 0.0;
                vvxy[ix][iy] = 0.0;
                vvyx[ix][iy] = 0.0;
                vvyy[ix][iy] = 0.0;

                #pragma omp parallel for
                for (int jx = 0; jx <= Lx; jx++) {
                    int kx = ix - jx + Nx;
                    if (kx < 0) continue;
                    if (kx > Lx) continue;
                    for (int jy = 0; jy <= Ly; jy++) {
                        int ky = iy - jy + Ny;
                        if (ky < 0) continue;
                        if (ky > Ly) continue;

                        if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny)) == 0) {
                            continue;
                        }

                        //vv
                        vvxx[ix][iy] += vx[jx][jy] * vx[kx][ky] * DVQ;
                        vvxy[ix][iy] += vx[jx][jy] * vy[kx][ky] * DVQ;
                        vvyx[ix][iy] += vy[jx][jy] * vx[kx][ky] * DVQ;
                        vvyy[ix][iy] += vy[jx][jy] * vy[kx][ky] * DVQ;
                    }
                }
            }
        }


        for (int ix = 0; ix <= Lx; ix++) {
            for (int iy = 0; iy <= Ly; iy++) {

                double Qx = -Qx0 + ix * DQx;
                double Qy = -Qy0 + iy * DQy;

                double Q2 = Qx * Qx + Qy * Qy;
                if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny)) == 0) {
                    continue;
                }

                double Q = sqrt(Q2);
                double nx = Qx / Q;
                double ny = Qy / Q;

                double fnvx =
                        (vvxx[ix][iy]) * Qx +
                        (vvxy[ix][iy]) * Qy;
                double fnvy =
                        (vvyx[ix][iy]) * Qx +
                        (vvyy[ix][iy]) * Qy;


                double flvx = -nu * Q2 * vx[ix][iy];
                double flvy = -nu * Q2 * vy[ix][iy];


                double fvx = flvx + fnvx;
                double fvy = flvy + fnvy;

                double fvn = fvx * nx + fvy * ny;

                fvx -= fvn * nx;
                fvy -= fvn * ny;

                vx[ix][iy] += fvx * dt;
                vy[ix][iy] += fvy * dt;

            }

            double EE = 0.0;
            for (int ix = 0; ix <= Lx; ix++) {
                for (int iy = 0; iy <= Ly; iy++) {
                    EE += vx[ix][iy] * vx[ix][iy] + vy[ix][iy] * vy[ix][iy];
                }
            }
            printf("%i %g\n", it, EE);

        }
    }


    return 0;
}
