#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>


const double vc = 0.1;
const double vphi = 0.2;
const double E0= 0.045;

int main() {
    double nu = 0.000001;

    double Elast[100];
    int pos=0;

    double Qx0 = 3.0;
    double Qy0 = 3.0;

    int Nx = 8;
    int Ny = 8;

    int Nt = 10000000;
//    double dt = 0.0025;
    double dt = 0.01;

    int NDimx = 2 * Nx + 4, NDimy = 2 * Ny + 4;

    double vx[NDimx][NDimy];
    double vy[NDimx][NDimy];
    double vz[NDimx][NDimy];

    int Lx = 2 * Nx;
    int Ly = 2 * Ny;

    double DQx = Qx0 / Nx;
    double DQy = Qy0 / Ny;

    double DVQ = DQx * DQy / ((2 * M_PI) * (2 * M_PI));

    double vvxx[NDimx][NDimy], vvxy[NDimx][NDimy];
    double vvyx[NDimx][NDimy], vvyy[NDimx][NDimy];

    for (int ix = 0; ix <= Lx; ix++) {
        for (int iy = 0; iy <= Ly; iy++) {
            double Qx = -Qx0 + ix * DQx;
            double Qy = -Qy0 + iy * DQy;

            double Q2 = Qx * Qx + Qy * Qy;
            double norm = 1 / sqrt(Q2);


            if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny)) == 0) {
                continue;
            }

            vx[ix][iy] = vphi * -Qy * norm;
            vy[ix][iy] = vphi * Qx * norm;
        }
    }

    vx[Nx + 1][Ny] = 0.0;
    vy[Nx + 1][Ny] = vc;


    vx[Nx][Ny + 1] = -vc;
    vy[Nx][Ny + 1] = 0.0;


    vx[Nx - 1][Ny] = 0.0;
    vy[Nx - 1][Ny] = -vc;

    vx[Nx][Ny - 1] = vc;
    vy[Nx][Ny - 1] = 0.0;



    for (int it = 1; it <= Nt; it++) {

	#pragma omp for
        for (int ix = 0; ix <= Lx; ix++) {
            for (int iy = 0; iy <= Ly; iy++) {

                vvxx[ix][iy] = 0.0;
                vvxy[ix][iy] = 0.0;
                vvyx[ix][iy] = 0.0;
                vvyy[ix][iy] = 0.0;


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

                if (abs(ix - Nx) == 1 && iy == Ny) {
                    continue;
                }

                if (abs(iy - Ny) == 1 && ix == Nx) {
                    continue;
                }

                double Q = sqrt(Q2);
                double nx = Qx / Q;
                double ny = Qy / Q;

                double fnvx = vvxx[ix][iy] * Qx + vvxy[ix][iy] * Qy;
                double fnvy = vvyx[ix][iy] * Qx + vvyy[ix][iy] * Qy;

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
        }

        double EE = 0.0;
        double vphi = 0.0;
        for (int ix = 0; ix <= Lx; ix++) {
            for (int iy = 0; iy <= Ly; iy++) {
                double Qx = -Qx0 + ix * DQx;
                double Qy = -Qy0 + iy * DQy;

                if (abs(ix - Nx) == 1 && iy == Ny) {
                    continue;
                }

                if (abs(iy - Ny) == 1 && ix == Nx) {
                    continue;
                }

//                    double Qz = -Qz0 + iz * DQz;
                EE += (vx[ix][iy] * vx[ix][iy] + vy[ix][iy] * vy[ix][iy])* DVQ;
                double vphi = (-Qy * vx[ix][iy] + Qx * vy[ix][iy]) / sqrt(Qx * Qx + Qy * Qy);
                printf("%f %f %f\n", Qx, Qy, vphi);

            }
        }
	if ( it < 100 ) 
	{
	nu = nu*(EE/E0);
	}
	else if ( it > 100 && it < 650000 )
	{
	Elast[pos]=EE;
        pos = pos + 1;
        if (pos == 10)
		{
                double Eav=0;
		for (int i=0;i<=10;++i)
		    {
		    Eav=Eav+Elast[i];
		    }
		Eav=Eav/10;		
		nu = nu*(Eav/E0);
		pos=0;
		}
	}
        else {
	nu = nu*(EE/E0);	
	}
puts("============================");
	
//        printf("%i %g\n", it, EE);
    }

    return 0;
}
