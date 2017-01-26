#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{

	double Qx0=6.0, Qy0=6.0, Qz0=6.0;
	int  Nx=1, Ny=1, Nz=10;
	int NDimx=2*Nx+2, NDimy=2*Ny+2, NDimz=2*Nz+2;
	double vx[NDimx][NDimy][NDimz], vy[NDimx][NDimy][NDimz], vz[NDimx][NDimy][NDimz],
	bx[NDimx][NDimy][NDimz], by[NDimx][NDimy][NDimz], bz[NDimx][NDimy][NDimz],
	DV[NDimx][NDimy][NDimz];
	int Lx=2*Nx, Ly= 2*Ny, Lz=2*Nz;
	double DQx=Qx0/Nx, DQy=Qy0/Ny, DQz=Qz0/Nz, rnx, rny, rnz, vn, bn;
	double DVQ=DQx*DQy*DQz/((2*M_PI)*(2*M_PI)*(2*M_PI));
	double alpha=0.0;
	double betax,betaz;
	double betay=0.0;
	double omega=-2./3.;
	double fnvx,fnvy,fnvz,fnbx,fnby,fnbz;
	double Q2,Qbeta;

	double fvx[NDimx][NDimy][NDimz];
	double fvy[NDimx][NDimy][NDimz];
	double fvz[NDimx][NDimy][NDimz];
	double fbx[NDimx][NDimy][NDimz];
	double fby[NDimx][NDimy][NDimz];
	double fbz[NDimx][NDimy][NDimz];
	double rnuk=0.0;
	double rnum=0.0;
	double Tempi,Nt,dt,t;
	int it;


  double Ei[NDimx][NDimy][NDimz];
	double Ef[NDimx][NDimy][NDimz];

	betax = sin(alpha);
	betaz = cos(alpha);

    double Q=0.0;
// Inicializacija
	for (int ix=0;ix<=Lx; ix++)
	{double Qx=-Qx0+ix*DQx;
		for (int iy=0;iy<=Ly; iy++)
		{double Qy=-Qy0+iy*DQy;
			for (int iz=0;iz<=Lz; iz++)
			{double Qz=-Qz0+iz*DQz;

				vx[ix][iy][iz] = 1.0;
				vy[ix][iy][iz] = 1.0;
				vz[ix][iy][iz] = 0.0;

				bx[ix][iy][iz] = 1.0;
				by[ix][iy][iz] = 1.0;
				bz[ix][iy][iz] = 0.0;

				 Q2=Qx*Qx+Qy*Qy+Qz*Qz;
				 if ( ((ix - Nx)*(ix - Nx) + (iy-Ny)*(iy-Ny) + (iz-Nz)*(iz-Nz)) == 0 )
				 {
					 continue;
				 }
				 Q=sqrt(Q2);
				 rnx=Qx/Q; rny=Qy/Q; rnz=Qz/Q;

				 vn = vx[ix][iy][iz]*rnx + vy[ix][iy][iz]*rny + vz[ix][iy][iz]*rnz;
				 bn = bx[ix][iy][iz]*rnx + by[ix][iy][iz]*rny + bz[ix][iy][iz]*rnz;

				vx[ix][iy][iz]-=vn*rnx;
				vy[ix][iy][iz]-=vn*rny;
				vz[ix][iy][iz]-=vn*rnz;

				bx[ix][iy][iz]-=bn*rnx;
				by[ix][iy][iz]-=bn*rny;
				bz[ix][iy][iz]-=bn*rnz;

				// printf("%g %g %g %g %g %g\n",vvxx[ix][iy][iz],bbxx[ix][iy][iz],vvxy[ix][iy][iz],bbxy[ix][iy][iz],vvxz[ix][iy][iz],bbxz[ix][iy][iz]);
				Ei[ix][iy][iz] = vx[ix][iy][iz]*vx[ix][iy][iz] + vy[ix][iy][iz]*vy[ix][iy][iz]  + bx[ix][iy][iz]*bx[ix][iy][iz] + by[ix][iy][iz]*by[ix][iy][iz];

			}
		}
	}

		Tempi=25.0; Nt=1000; dt=Tempi/Nt;
        // nachalo na integrirane po vreme; tuk se vyrti
		for(it=1;it<=Nt; it++)
		{t=it*dt;
			for(int ix=0;ix<=Lx; ix++)
			{double Qx=-Qx0+ix*DQx;
				for(int iy=0;iy<=Ly; iy++)
				{double Qy=-Qy0+iy*DQy;
					for(int iz=0;iz<=Lz; iz++)
					{double Qz=-Qz0+iz*DQz;

					 double fvn,fbn,aaa,fnvx,fnvy,fnvz,fnbx,fnby,fnbz,flvx,flvy,flvz,flbx,flby,flbz;

					 Q2=Qx*Qx+Qy*Qy+Qz*Qz;
					 if ( ((ix - Nx)*(ix - Nx) + (iy-Ny)*(iy-Ny) + (iz-Nz)*(iz-Nz)) == 0 )
					 {
						 continue;
					 }
					 Q=sqrt(Q2);
					 rnx=Qx/Q; rny=Qy/Q; rnz=Qz/Q;
					 Qbeta=Qx*betax+Qz*betaz; // betay=0.



					 aaa = 2.0*( - omega*(rnx*vy[ix][iy][iz]-rny*vx[ix][iy][iz]));
					 flvx =                  aaa*rnx - 2.0*omega*vy[ix][iy][iz] + Qbeta*bx[ix][iy][iz] - rnuk*Q2*vx[ix][iy][iz];
					 flvy = aaa*rny + 2.0*omega*vx[ix][iy][iz] + Qbeta*by[ix][iy][iz] - rnuk*Q2*vy[ix][iy][iz];
					 flvz =                  aaa*rnz                            + Qbeta*bz[ix][iy][iz] - rnuk*Q2*vz[ix][iy][iz];

					 flbx =      	       - Qbeta*vx[ix][iy][iz] - rnum*Q2*vx[ix][iy][iz];
					 flby = bx[ix][iy][iz] - Qbeta*vy[ix][iy][iz] - rnum*Q2*vy[ix][iy][iz];
					 flbz =      	       - Qbeta*vz[ix][iy][iz] - rnum*Q2*vz[ix][iy][iz];


					 fvn = flvx*rnx + flvy*rny + flvz*rnz;
					 fbn = flbx*rnx + flby*rny + flbz*rnz;
					 printf("Qx=%g \t Qy=%g \t Qz=%g \t fvn=%g \t fbn=%g \n",Qx,Qy,Qz,fvn,fbn);

					 flvx-= fvn*rnx;
					 flvy-= fvn*rny;
					 flvz-= fvn*rnz;

					 flbx-= fbn*rnx;
					 flby-= fbn*rny;
					 flbz-= fbn*rnz;

					 vx[ix][iy][iz]+=flvx*dt; vy[ix][iy][iz]+=flvy*dt; vz[ix][iy][iz]+=flvz*dt;
           				 bx[ix][iy][iz]+=flbx*dt; by[ix][iy][iz]+=flby*dt; bz[ix][iy][iz]+=flbz*dt;


					Ef[ix][iy][iz] = vx[ix][iy][iz]*vx[ix][iy][iz] + vy[ix][iy][iz]*vy[ix][iy][iz] + bx[ix][iy][iz]*bx[ix][iy][iz] + by[ix][iy][iz]*by[ix][iy][iz];
					}//iz
				}//iy
			}//ix

		}//it



							for(int iz=0;iz<=Lz; iz++)
							{
							 double Qz=-Qz0+iz*DQz;
							 printf("%g %g\n",Qz, pow( (1/(2.0*Tempi))*log(Ef[1][1][iz]/Ei[1][1][iz]),2) );
							}





	return 0;
}// main
