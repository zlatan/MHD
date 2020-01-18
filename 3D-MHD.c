#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

#define sign(z) ( ( (z) < 0 )  ?  -1   : ( (z) > 0 ) )

const double nu = 0.0022;
const double vc = 2.5;
const double vphi = 0.2;
const double alpha = 0.0;
const double nu_k = 0.001;
const double nu_m = 0.001;
const double omega = -2. / 3.;

int main() {
    double Qx0 = 4.0;
    double Qy0 = 4.0;
    double Qz0 = 4.0;

    int Nx = 4;
    int Ny = 4;
    int Nz = 4;

    int Nt = 100000;
    double dt = 0.0025;
    
    int NDimx = 2 * Nx + 4, NDimy = 2 * Ny + 4, NDimz = 2 * Nz + 4;
    
    double vx[NDimx][NDimy][NDimz];
    double vy[NDimx][NDimy][NDimz];
    double vz[NDimx][NDimy][NDimz];
    
    double bx[NDimx][NDimy][NDimz];
    double by[NDimx][NDimy][NDimz];
    double bz[NDimx][NDimy][NDimz];

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

                bx[ix][iy][iz] = vx[ix][iy][iz]*;
                by[ix][iy][iz] = 0.0;
                bz[ix][iy][iz] = 0.0;
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

                    double fnvx =(vvxx[ix][iy][iz] + bbxx[ix][iy][iz]) * Qx + (vvxy[ix][iy][iz] + bbxy[ix][iy][iz]) * Qy +(vvxz[ix][iy][iz] + bbxz[ix][iy][iz]) * Qz;
                    double fnvy =(vvyx[ix][iy][iz] + bbyx[ix][iy][iz]) * Qx + (vvyy[ix][iy][iz] + bbyy[ix][iy][iz]) * Qy +(vvyz[ix][iy][iz] + bbyz[ix][iy][iz]) * Qz;
                    double fnvz =(vvzx[ix][iy][iz] + bbzx[ix][iy][iz]) * Qx + (vvzy[ix][iy][iz] + bbzy[ix][iy][iz]) * Qy +(vvzz[ix][iy][iz] + bbzz[ix][iy][iz]) * Qz;

                    double fnbx =(bvxx[ix][iy][iz] - vbxx[ix][iy][iz]) * Qx + (bvxy[ix][iy][iz] - vbxy[ix][iy][iz]) * Qy +(bvxz[ix][iy][iz] - vbxz[ix][iy][iz]) * Qz;
                    double fnby =(bvyx[ix][iy][iz] - vbyx[ix][iy][iz]) * Qx + (bvyy[ix][iy][iz] - vbyy[ix][iy][iz]) * Qy +(bvyz[ix][iy][iz] - vbyz[ix][iy][iz]) * Qz;
                    double fnbz =(bvzx[ix][iy][iz] - vbzx[ix][iy][iz]) * Qx + (bvzy[ix][iy][iz] - vbzy[ix][iy][iz]) * Qy +(bvzz[ix][iy][iz] - vbzz[ix][iy][iz]) * Qz;

                    double ct = 2.0 * (ny * vx[ix][iy][iz] - omega * (nx * vy[ix][iy][iz] - ny * vx[ix][iy][iz]));
                    double flvx = ct * nx + 2.0 * omega * vy[ix][iy][iz] + Qalpha * bx[ix][iy][iz] -
                                  nu_k * Q2 * vx[ix][iy][iz];
                    double flvy = -vx[ix][iy][iz] + ct * ny - 2.0 * omega * vx[ix][iy][iz] + Qalpha * by[ix][iy][iz] -
                                  nu_k * Q2 * vy[ix][iy][iz];
                    double flvz = ct * nz + Qalpha * bz[ix][iy][iz] - nu_k * Q2 * vz[ix][iy][iz];

                    double flbx = -Qalpha * vx[ix][iy][iz] - nu_m * Q2 * bx[ix][iy][iz];
                    double flby = bx[ix][iy][iz] - Qalpha * vy[ix][iy][iz] - nu_m * Q2 * by[ix][iy][iz];
                    double flbz = -Qalpha * vz[ix][iy][iz] - nu_m * Q2 * bz[ix][iy][iz];

                    double fvx = flvx + fnvx;
                    double fvy = flvy + fnvy;
                    double fvz = flvz + fnvz;

                    double fbx = flbx + fnbx;
                    double fby = flby + fnby;
                    double fbz = flbz + fnbz;

                    double fvn = fvx * nx + fvy * ny + fvz * nz;
                    double fbn = fbx * nx + fby * ny + fbz * nz;

                    fvx -= fvn * nx;
                    fvy -= fvn * ny;
                    fvz -= fvn * nz;
                    fbx -= fbn * nx;
                    fby -= fbn * ny;
                    fbz -= fbn * nz;

                    vx[ix][iy][iz] += fvx * dt;
                    vy[ix][iy][iz] += fvy * dt;
                    vz[ix][iy][iz] += fvz * dt;

                    bx[ix][iy][iz] += fbx * dt;
                    by[ix][iy][iz] += fby * dt;
                    bz[ix][iy][iz] += fbz * dt;

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
                    EE += (vx[ix][iy][iz] * vx[ix][iy][iz] + vy[ix][iy][iz] * vy[ix][iy][iz] +
                          vz[ix][iy][iz] * vz[ix][iy][iz])*DVQ;

//                double vphi = (-Qy * vx[ix][iy][Nz] + Qx * vy[ix][iy][Nz]) / sqrt(Qx * Qx + Qy * Qy);
//                printf("%f %f %f\n", Qx, Qy, vphi);
                }
            }
        }
            printf("%i %.19f\n", it, EE);

    }

    return 0;
}