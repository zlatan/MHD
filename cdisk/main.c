#include <stdio.h>
#include <math.h>

const double alpha = 0.0;
const double omega = -2. / 3.;
const double rnuk = 0.0;
const double rnum = 0.0;


int main() {


    double Qx0 = 2.0;
    double Qy0 = 2.0;
    double Qz0 = 2.0;

    int Nx = 4;
    int Ny = 4;
    int Nz = 4;

    double Tempi = 25.0;
    int Nt = 1000;
    double dt = Tempi / Nt;

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

    double Ei[NDimx][NDimy][NDimz];
    double Ef[NDimx][NDimy][NDimz];

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


    FILE *f = fopen("q_v.txt", "w");

    for (int ix = 0; ix <= Lx; ix++) {
        for (int iy = 0; iy <= Ly; iy++) {
            for (int iz = 0; iz <= Lz; iz++) {
                double Qx = -Qx0 + ix * DQx;
                double Qy = -Qy0 + iy * DQy;
                double Qz = -Qz0 + iz * DQz;

                vx[ix][iy][iz] = 1.0;
                vy[ix][iy][iz] = 1.0;
                vz[ix][iy][iz] = 0.0;

                bx[ix][iy][iz] = 1.0;
                by[ix][iy][iz] = 1.0;
                bz[ix][iy][iz] = 0.0;


                double Q2 = Qx * Qx + Qy * Qy + Qz * Qz;
                double Q = sqrt(Q2);

                if (((ix - Nx) * (ix - Nx) + (iy - Ny) * (iy - Ny) + (iz - Nz) * (iz - Nz)) == 0) {
                    continue;
                }
                double nx = Qx / Q;
                double ny = Qy / Q;
                double nz = Qz / Q;

                double vn = vx[ix][iy][iz] * nx + vy[ix][iy][iz] * ny + vz[ix][iy][iz] * nz;
                double bn = bx[ix][iy][iz] * nx + by[ix][iy][iz] * ny + bz[ix][iy][iz] * nz;

                vx[ix][iy][iz] -= vn * nx;
                vy[ix][iy][iz] -= vn * ny;
                vz[ix][iy][iz] -= vn * nz;

                bx[ix][iy][iz] -= bn * nx;
                by[ix][iy][iz] -= bn * ny;
                bz[ix][iy][iz] -= bn * nz;

                bbxx[ix][iy][iz] = 0.0;
                bbxy[ix][iy][iz] = 0.0;
                bbxz[ix][iy][iz] = 0.0;
                bbyx[ix][iy][iz] = 0.0;
                bbyy[ix][iy][iz] = 0.0;
                bbyz[ix][iy][iz] = 0.0;
                bbzx[ix][iy][iz] = 0.0;
                bbzy[ix][iy][iz] = 0.0;
                bbzz[ix][iy][iz] = 0.0;
                //vv
                vvxx[ix][iy][iz] = 0.0;
                vvxy[ix][iy][iz] = 0.0;
                vvxz[ix][iy][iz] = 0.0;
                vvyx[ix][iy][iz] = 0.0;
                vvyy[ix][iy][iz] = 0.0;
                vvyz[ix][iy][iz] = 0.0;
                vvzx[ix][iy][iz] = 0.0;
                vvzy[ix][iy][iz] = 0.0;
                vvzz[ix][iy][iz] = 0.0;
                //bv
                bvxx[ix][iy][iz] = 0.0;
                bvxy[ix][iy][iz] = 0.0;
                bvxz[ix][iy][iz] = 0.0;
                bvyx[ix][iy][iz] = 0.0;
                bvyy[ix][iy][iz] = 0.0;
                bvyz[ix][iy][iz] = 0.0;
                bvzx[ix][iy][iz] = 0.0;
                bvzy[ix][iy][iz] = 0.0;
                bvzz[ix][iy][iz] = 0.0;
                //vb
                vbxx[ix][iy][iz] = 0.0;
                vbxy[ix][iy][iz] = 0.0;
                vbxz[ix][iy][iz] = 0.0;
                vbyx[ix][iy][iz] = 0.0;
                vbyy[ix][iy][iz] = 0.0;
                vbyz[ix][iy][iz] = 0.0;
                vbzx[ix][iy][iz] = 0.0;
                vbzy[ix][iy][iz] = 0.0;
                vbzz[ix][iy][iz] = 0.0;
            }
        }
    }


    for (int it = 1; it <= Nt; it++) {
        puts("debug");
        for (int ix = 0; ix <= Lx; ix++) {
            for (int iy = 0; iy <= Ly; iy++) {
                for (int iz = 0; iz <= Lz; iz++) {
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

                                double Qx = -Qx0 + ix * DQx;
                                double Qy = -Qy0 + iy * DQy;
                                double Qz = -Qz0 + iz * DQz;

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

//                                    fprintf(f, " %+d %+d %+d %+d %+d %+d %+.8f\n", ix, iy, iz, jx, jy , jz, bbxx[ix][iy][iz]);
//                                    fprintf(f, " %+d %+d %+d %+d %+d %+d %+.8f %+.8f %+.8f\n", ix, iy, iz, jx, jy , jz, bbxx[ix][iy][iz], bbxy[ix][iy][iz],bbxz[ix][iy][iz]);

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

                    double Q = sqrt(Q2);
                    double nx = Qx / Q;
                    double ny = Qy / Q;
                    double nz = Qz / Q;
                    double Qalpha = Qy * sin(alpha) + Qz * cos(alpha);

                    double fnvx = (vvxx[ix][iy][iz] + bbxx[ix][iy][iz]) * Qx +  (vvxy[ix][iy][iz] + bbxy[ix][iy][iz]) * Qy + (vvxz[ix][iy][iz] + bbxz[ix][iy][iz]) * Qz;
                    double fnvy =
                            (vvyx[ix][iy][iz] + bbyx[ix][iy][iz]) * Qx +
                            (vvyy[ix][iy][iz] + bbyy[ix][iy][iz]) * Qy +
                            (vvyz[ix][iy][iz] + bbyz[ix][iy][iz]) * Qz;
                    double fnvz =
                            (vvzx[ix][iy][iz] + bbzx[ix][iy][iz]) * Qx +
                            (vvzy[ix][iy][iz] + bbzy[ix][iy][iz]) * Qy +
                            (vvzz[ix][iy][iz] + bbzz[ix][iy][iz]) * Qz;


                    double fnbx =
                            (bvxx[ix][iy][iz] - vbxx[ix][iy][iz]) * Qx +
                            (bvxy[ix][iy][iz] - vbxy[ix][iy][iz]) * Qy +
                            (bvxz[ix][iy][iz] - vbxz[ix][iy][iz]) * Qz;
                    double fnby =
                            (bvyx[ix][iy][iz] - vbyx[ix][iy][iz]) * Qx +
                            (bvyy[ix][iy][iz] - vbyy[ix][iy][iz]) * Qy +
                            (bvyz[ix][iy][iz] - vbyz[ix][iy][iz]) * Qz;
                    double fnbz =
                            (bvzx[ix][iy][iz] - vbzx[ix][iy][iz]) * Qx +
                            (bvzy[ix][iy][iz] - vbzy[ix][iy][iz]) * Qy +
                            (bvzz[ix][iy][iz] - vbzz[ix][iy][iz]) * Qz;


                    double ct = 2.0 * (ny * vx[ix][iy][iz] - omega * (nx * vy[ix][iy][iz] - ny * vx[ix][iy][iz]));
                    double flvx = ct * nx + 2.0 * omega * vy[ix][iy][iz] + Qalpha * bx[ix][iy][iz] -
                                  rnuk * Q2 * vx[ix][iy][iz];
                    double flvy = -vx[ix][iy][iz] + ct * ny - 2.0 * omega * vx[ix][iy][iz] + Qalpha * by[ix][iy][iz] -
                                  rnuk * Q2 * vy[ix][iy][iz];
                    double flvz = ct * nz + Qalpha * bz[ix][iy][iz] - rnuk * Q2 * vz[ix][iy][iz];

                    double flbx = -Qalpha * vx[ix][iy][iz] - rnum * Q2 * vx[ix][iy][iz];
                    double flby = bx[ix][iy][iz] - Qalpha * vy[ix][iy][iz] - rnum * Q2 * vy[ix][iy][iz];
                    double flbz = -Qalpha * vz[ix][iy][iz] - rnum * Q2 * vz[ix][iy][iz];

                        double fvx = flvx + fnvx;
                        double fvy = flvy + fnvy;
                        double fvz = flvz + fnvz;

                        double fbx = flbx + fnbx;
                        double fby = flby + fnby;
                        double fbz = flbz + fnbz;

//                    fprintf(f, " %+.1f %+.1f %+.1f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f %+.8f\n", Qx, Qy, Qz, bbxx[ix][iy][iz], bbxy[ix][iy][iz], bbxz[ix][iy][iz], bbyx[ix][iy][iz], bbyy[ix][iy][iz], bbyz[ix][iy][iz], bbzx[ix][iy][iz], bbzy[ix][iy][iz], bbzz[ix][iy][iz]);
                    fprintf(f, " %+.1f %+.1f %+.1f %+.16f\n", Qx, Qy, Qz, vvxx[ix][iy][iz]);


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

//                        fprintf(f, " %+.1f %+.1f %+.1f %+.2f %+.2f %+.2f\n", Qx, Qy, Qz, vx[ix][iy][iz], vy[ix][iy][iz],vz[ix][iy][iz]);
                }
            }
        }
        break;
    }
    fclose(f);
    printf("Hello, World!\n");


    return 0;
}