#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{
	double A=1.0;
	double C=1.0;

	double Qx0=3.0, Qy0=3.0, Qz0=3.0;
	int  Nz=100;
	int  NDim=2*Nz+2;
	double vx[NDim], vy[NDim], vz[NDim];
	double bx[NDim], by[NDim], bz[NDim];
	double bbzz[NDim], vvzz[NDim], bvzz[NDim], vbzz[NDim];
	double v[NDim];
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

	double f;
	double fvx[NDim];
	double fvy[NDim];
	double fvz[NDim];
	double fbx[NDim];
	double fby[NDim];
	double fbz[NDim];
	double rnuk=0.1;
	double rnum=0.1;
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

		Tempi=25.0; Nt=20; dt=Tempi/Nt;
//	 Tempi=25.0; Nt=100; dt=Tempi/Nt;

        // nachalo na integrirane po vreme; tuk se vyrti
		for(it=1;it<=Nt; it++)
		{t=it*dt;


////////////////////////////////////////////////////////////////////

						for (int iz=0;iz<=Lz; iz++)
								{double Qz=-Qz0+iz*DQ;
									// summation cycles for nonlinear terms
											for (int jz=0;jz<=Lz; jz++)
											{// razgele, naj setne zapoichvame rabota
											 int kz=iz-jz+Nz; if(kz<0) continue; if(kz>Lz) continue;
											 bbzz[iz]+=bz[jz]*bz[kz]*DV[jz];
											 vvzz[iz]+=vz[jz]*vz[kz]*DV[jz];
											 bvzz[iz]+=bz[jz]*vz[kz]*DV[jz];
										   vbzz[iz]+=vz[jz]*bz[kz]*DV[jz];
											}//kz
								}//iz

///////////////////////////////////////

					for(int iz=0;iz<=Lz; iz++)
					{double Qz=-Qz0+iz*DQ;

					 double fvn,fbn,aaa,fnvx,fnvy,fnvz,fnbx,fnby,fnbz,flvx,flvy,flvz,flbx,flby,flbz;


					 Qbeta=Qz*betaz; // betay=0.


					 fnvz = (vvzz[iz] + bbzz[iz])*Qz;
					 fnbz = (bvzz[iz] - vbzz[iz])*Qz;

					 fvn = fnvz;
					 fbn = fnbz;

					 fnvz-= fvn;
  				 fnbz-= fbn;


					 flvx =           + 2.0*omega*vy[iz] + Qz*bx[iz] - rnuk*Qz*Qz*vx[iz];
					 flvy =-vx[iz]    - 2.0*omega*vx[iz] + Qz*by[iz] - rnuk*Qz*Qz*vy[iz];
					 flvz =             		               Qz*bz[iz] - rnuk*Qz*Qz*vz[iz];

					 flbx =        - Qz*vx[iz] - rnum*Qz*Qz*bx[iz];
					 flby = bx[iz] - Qz*vy[iz] - rnum*Qz*Qz*by[iz];
					 flbz =        - Qz*vz[iz] - rnum*Qz*Qz*bz[iz];


///////////////////////////////////////////////////////////////////////////////////////

flvz-= fvn;
flbz-= fbn;

////////////////////////////////////////////////////////////////////////////////////////


					 vx[iz]+=flvx*dt;					 vy[iz]+=flvy*dt;					 vz[iz]+=(flvz+fnvz)*dt;
			     bx[iz]+=flbx*dt;					 by[iz]+=flby*dt;					 bz[iz]+=(flbz+fnbz)*dt;

					 Ef[iz] = vx[iz]*vx[iz] + vy[iz]*vy[iz] + bx[iz]*bx[iz] + by[iz]*by[iz];
					}//iz
				 //printf("%g %g\n",t,bx[0]);

		}//it


							for(int iz=0;iz<=Lz; iz++)
							{
							 double Qz=-Qz0+iz*DQ;
							//printf("%g %g\n",Qz,v[iz]);
							//  printf("%g %g\n",Qz, (1/(2.0*Tempi))*log(Ef[iz]/Ei[iz]) );
							printf("%g %g\n",Qz,Ef[iz]);
							}

	//printf("%g\n",(-1/Tempi)*log(v[Lz]/v[0]));
// p "../../vasko/MHD/disk3.out" w l

	return 0;
}// main
