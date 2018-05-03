#include <stdio.h>
#include <math.h>

const double nu_k=0.001;
const double nu_m=0.001;
const double omega=-2.0/3.0;
const double alpha=0.0;

typedef struct Vector
{
    double x,y,z;
} Vector;


typedef struct TriVector
{
    Vector wave_vector;
    Vector velocity;
    Vector magnetic_field;

} TriVector;

typedef struct Energy
{
    Vector wave_vector;
    double energy;
}Energy;


TriVector local_linear(TriVector mhd)
{
    double Q2 = mhd.wave_vector.x*mhd.wave_vector.x+ mhd.wave_vector.y*mhd.wave_vector.y+mhd.wave_vector.z*mhd.wave_vector.z;
    double Q = sqrt(Q2);
    double Qa = mhd.wave_vector.y*sin(alpha) + mhd.wave_vector.z*cos(alpha);
    double nx = mhd.wave_vector.x/Q;
    double ny = mhd.wave_vector.y/Q;
    double nz = mhd.wave_vector.z/Q;

    TriVector next_mhd;

    double cm = 2.0*(ny*mhd.velocity.x - omega*(nx*mhd.velocity.y-ny*mhd.velocity.x));
    next_mhd.velocity.x =                  cm*nx + 2.0*omega*mhd.velocity.y + Qa*mhd.magnetic_field.x - nu_k*Q2*mhd.velocity.x;
    next_mhd.velocity.y =-mhd.velocity.x + cm*ny - 2.0*omega*mhd.velocity.x + Qa*mhd.magnetic_field.y - nu_k*Q2*mhd.velocity.y;
    next_mhd.velocity.z =                  cm*nz                            + Qa*mhd.magnetic_field.z - nu_k*Q2*mhd.velocity.z;

    next_mhd.magnetic_field.x =                      - Qa*mhd.velocity.x - nu_m*Q2*mhd.magnetic_field.x;
    next_mhd.magnetic_field.y = mhd.magnetic_field.x - Qa*mhd.velocity.y - nu_m*Q2*mhd.magnetic_field.y;
    next_mhd.magnetic_field.z =                      - Qa*mhd.velocity.z - nu_m*Q2*mhd.magnetic_field.z;

    next_mhd.wave_vector = mhd.wave_vector;

    return next_mhd;
}

Energy energy(TriVector mhd)
{
    Energy e;
    e.energy = mhd.velocity.x * mhd.velocity.x + mhd.velocity.y*mhd.velocity.y + mhd.velocity.z*mhd.velocity.z
             + mhd.magnetic_field.x*mhd.magnetic_field.x + mhd.magnetic_field.y*mhd.magnetic_field.y + mhd.magnetic_field.z*mhd.magnetic_field.z;
    e.wave_vector = mhd.wave_vector;
    return e;
}


int main() {
    TriVector init[1200];
    Energy init_e[1200];

    TriVector second[1200];
    Energy fin_e[1200];

    int c = 0;

    for (double x = 0.01; x <= 1.0; x = x + 0.1)
        for (double y = 0.01; y <= 1.0; y = y + 0.1)
            for (double z = 0.01; z <= 1.0; z = z + 0.1) {
                Vector Q, B, V;
                TriVector mhd;
                Q = (Vector) {x, y, z};
                B = (Vector) {1.0, 1.0, 0.0};
                V = (Vector) {1.0, 1.0, 0.0};
                mhd = (TriVector) {Q, V, B};
                init[c] = mhd;
                init_e[c] = energy(mhd);
                c++;
            }
//    local_linear(mhd);

    FILE *f = fopen("file.txt", "w");

    for (int i = 0; i <= 1000; i++)
    {
        second[i] = local_linear(init[i]);
        fin_e[i] = energy(second[i]);
        printf("%g %g %g \n",fin_e[i].wave_vector.y,fin_e[i].wave_vector.z,pow( (1/(2.0))*log(fin_e[i].energy/init_e[i].energy),2));
        fprintf(f, "%g %g %g \n",fin_e[i].wave_vector.y,fin_e[i].wave_vector.z,pow( (1/(2.0))*log(fin_e[i].energy/init_e[i].energy),2));

    }
    fclose(f);


    printf("Hello, World!\n");
    return 0;
}