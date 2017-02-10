#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{

	double Qx0=2.0, Qy0=2.0, Qz0=2.0;
	int  Nz=100;
	int  NDim=2*Nz+2;
	double vx[NDim], vy[NDim], vz[NDim];
	double bx[NDim], by[NDim], bz[NDim];
	double DV[NDim];
	int  Lz=2*Nz;
	double DQ=Qz0/Nz, rnx, rny, rnz, vn, bn;
	double DVQ=DQ/(2*M_PI);
	double alpha=0.0;
	double betax,betaz;
	double betay=0.0;
	double omega=-2./3.;
	double fnvx,fnvy,fnvz,fnbx,fnby,fnbz;
	double Q2,Qbeta;

	double fvx[NDim];
	double fvy[NDim];
	double fvz[NDim];
	double fbx[NDim];
	double fby[NDim];
	double fbz[NDim];
	double rnuk=0.0;
	double rnum=0.0;
	double Tempi,Nt,dt,t;
	int it;


  	double Ei[NDim];
	double Ef[NDim];

	betax = sin(alpha);
	betaz = cos(alpha);

    double Q=0.0;
// Inicializacija
			for (int iz=0;iz<=Lz; iz++)
			{
				double Qz=-Qz0+iz*DQ;

				vx[iz] = 1.0;
				vy[iz] = 1.0;
				vz[iz] = 0.0;

				bx[iz] = 1.0;
				by[iz] = 1.0;
				bz[iz] = 0.0;

				// printf("%g %g %g %g %g %g\n",vvxx[iz],bbxx[iz],vvxy[iz],bbxy[iz],vvxz[iz],bbxz[iz]);
				Ei[iz] = vx[iz]*vx[iz] + vy[iz]*vy[iz]  + bx[iz]*bx[iz] + by[iz]*by[iz];
				
			}

		Tempi=25.0; Nt=1000; dt=Tempi/Nt;
        // nachalo na integrirane po vreme; tuk se vyrti
		for(it=1;it<=Nt; it++)
		{t=it*dt;
					for(int iz=0;iz<=Lz; iz++)
					{double Qz=-Qz0+iz*DQ;

					 double fvn,fbn,aaa,fnvx,fnvy,fnvz,fnbx,fnby,fnbz,flvx,flvy,flvz,flbx,flby,flbz;


					 Qbeta=Qz*betaz; // betay=0.



					 flvx =         + 2.0*omega*vy[iz] + Qbeta*bx[iz];
					 flvy =-vx[iz]  - 2.0*omega*vx[iz] + Qbeta*by[iz];
					 flvz =             		   + Qbeta*bz[iz];

					 flbx =        - Qbeta*vx[iz];
					 flby = bx[iz] - Qbeta*vy[iz];
					 flbz =        - Qbeta*vz[iz];

					 vx[iz]+=flvx*dt; vy[iz]+=flvy*dt; vz[iz]+=flvz*dt;
           				 bx[iz]+=flbx*dt; by[iz]+=flby*dt; bz[iz]+=flbz*dt;


					Ef[iz] = vx[iz]*vx[iz] + vy[iz]*vy[iz] + bx[iz]*bx[iz] + by[iz]*by[iz];
					}//iz

		}//it



							for(int iz=0;iz<=Lz; iz++)
							{
							 double Qz=-Qz0+iz*DQ;
							 printf("%g %g\n",Qz, pow( (1/(2.0*Tempi))*log(Ef[iz]/Ei[iz]),2) );
							}





	return 0;
}// main
