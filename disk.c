#include <stdio.h>
#include <math.h>

int main(void)
{

	double Qx0=20., Qy0=3., Qz0=3.;
	Qx0=Qy0=Qz0=3.;// xxxx
	int Nx=30, Ny=6, Nz=6;
	Nx=10, Ny=10, Nz=10;//****
	int NDimx=2*Nx+2, NDimy=2*Ny+2, NDimz=2*Nz+2;
	double vx[NDimx][NDimy][NDimz], vy[NDimx][NDimy][NDimz], vz[NDimx][NDimy][NDimz],
	bx[NDimx][NDimy][NDimz], by[NDimx][NDimy][NDimz], bz[NDimx][NDimy][NDimz],
	DV[NDimx][NDimy][NDimz],
	bbyx[NDimx][NDimy][NDimz], bbyy[NDimx][NDimy][NDimz], bbyz[NDimx][NDimy][NDimz],
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
	vbzx[NDimx][NDimy][NDimz], vbzy[NDimx][NDimy][NDimz], vbzz[NDimx][NDimy][NDimz],
	rnvx[NDimx][NDimy][NDimz], rnvy[NDimx][NDimy][NDimz], rnvz[NDimx][NDimy][NDimz],
	rnbx[NDimx][NDimy][NDimz], rnby[NDimx][NDimy][NDimz], rnbz[NDimx][NDimy][NDimz];
	int Lx=2*Nx, Ly= 2*Ny, Lz=2*Nz;
	double DQx=Qx0/Nx, DQy=Qy0/Ny, DQz=Qz0/Nz, rnx, rny, rnz, vn, bn;
	double DVQ=DQx*DQy*DQz/((2*M_PI)*(2*M_PI)*(2*M_PI));
    double Q=0.0;

	for (int ix=0;ix<=Lx; ix++)
	{double Qx=-Qx0+ix*DQx;
		for (int iy=0;iy<=Ly; iy++)
		{double Qy=-Qy0+iy*DQy;
			for (int iz=0;iz<=Lz; iz++)
			{double Qz=-Qz0+iz*DQz;
		bbxx[ix][iy][iz]= bbxy[ix][iy][iz]= bbxz[ix][iy][iz]=
		bbyx[ix][iy][iz]= bbyy[ix][iy][iz]= bbyz[ix][iy][iz]=
		bbzx[ix][iy][iz]= bbzy[ix][iy][iz]= bbzz[ix][iy][iz]=
		vvxx[ix][iy][iz]= vvxy[ix][iy][iz]= vvxz[ix][iy][iz]=
		vvyx[ix][iy][iz]= vvyy[ix][iy][iz]= vvyz[ix][iy][iz]=
		vvzx[ix][iy][iz]= vvzy[ix][iy][iz]= vvzz[ix][iy][iz]=
		bvxx[ix][iy][iz]= bvxy[ix][iy][iz]= bvxz[ix][iy][iz]=
		bvyx[ix][iy][iz]= bvyy[ix][iy][iz]= bvyz[ix][iy][iz]=
		bvzx[ix][iy][iz]= bvzy[ix][iy][iz]= bvzz[ix][iy][iz]=
		vbxx[ix][iy][iz]= vbxy[ix][iy][iz]= vbxz[ix][iy][iz]=
		vbyx[ix][iy][iz]= vbyy[ix][iy][iz]= vbyz[ix][iy][iz]=
		vbzx[ix][iy][iz]= vbzy[ix][iy][iz]= vbzz[ix][iy][iz]=0.0;
		vx[ix][iy][iz]=exp(-(Qx*Qx+Qy*Qy+Qz*Qz));
		DV[ix][iy][iz]=DVQ;
			// summation cycles for nonlinear terms
				for (int jx=0;jx<=Lx; jx++)
				{int kx=ix-jx+Nx; if(kx<0) continue; if(kx>Lx) continue;
					for (int jy=0;jy<=Ly; jy++)
					{int ky=iy-jy+Ny; if(ky<0) continue; if(ky>Ly) continue;
						for (int jz=0;jz<=Lz; jz++)
						{// razgele, naj setne zapoichvame rabota
						 int kz=iz-jz+Nz; if(kz<0) continue; if(kz>Lz) continue;
						 //bb
							bbxx[ix][iy][iz]+=bx[jx][jy][jz]*bx[kx][ky][kz]*DV[jx][jy][jz]; bbxy[ix][iy][iz]+=bx[jx][jy][jz]*by[kx][ky][kz]*DV[jx][jy][jz]; bbxz[ix][iy][iz]+=bx[jx][jy][jz]*bz[kx][ky][kz]*DV[jx][jy][jz];
							bbyx[ix][iy][iz]+=by[jx][jy][jz]*bx[kx][ky][kz]*DV[jx][jy][jz]; bbyy[ix][iy][iz]+=by[jx][jy][jz]*by[kx][ky][kz]*DV[jx][jy][jz]; bbyz[ix][iy][iz]+=by[jx][jy][jz]*bz[kx][ky][kz]*DV[jx][jy][jz];
							bbzx[ix][iy][iz]+=bz[jx][jy][jz]*bx[kx][ky][kz]*DV[jx][jy][jz]; bbzy[ix][iy][iz]+=bz[jx][jy][jz]*by[kx][ky][kz]*DV[jx][jy][jz]; bbzz[ix][iy][iz]+=bz[jx][jy][jz]*bz[kx][ky][kz]*DV[jx][jy][jz];
						 //vv
							vvxx[ix][iy][iz]+=vx[jx][jy][jz]*vx[kx][ky][kz]*DV[jx][jy][jz]; vvxy[ix][iy][iz]+=vx[jx][jy][jz]*vy[kx][ky][kz]*DV[jx][jy][jz]; vvxz[ix][iy][iz]+=vx[jx][jy][jz]*vz[kx][ky][kz]*DV[jx][jy][jz];
							vvyx[ix][iy][iz]+=vy[jx][jy][jz]*vx[kx][ky][kz]*DV[jx][jy][jz]; vvyy[ix][iy][iz]+=vy[jx][jy][jz]*vy[kx][ky][kz]*DV[jx][jy][jz]; vvyz[ix][iy][iz]+=vy[jx][jy][jz]*vz[kx][ky][kz]*DV[jx][jy][jz];
							vvzx[ix][iy][iz]+=vz[jx][jy][jz]*vx[kx][ky][kz]*DV[jx][jy][jz]; vvzy[ix][iy][iz]+=vz[jx][jy][jz]*vy[kx][ky][kz]*DV[jx][jy][jz]; vvzz[ix][iy][iz]+=vz[jx][jy][jz]*vz[kx][ky][kz]*DV[jx][jy][jz];
						 //bv
							bvxx[ix][iy][iz]+=bx[jx][jy][jz]*vx[kx][ky][kz]*DV[jx][jy][jz]; bvxy[ix][iy][iz]+=bx[jx][jy][jz]*vy[kx][ky][kz]*DV[jx][jy][jz]; bvxz[ix][iy][iz]+=bx[jx][jy][jz]*vz[kx][ky][kz]*DV[jx][jy][jz];
							bvyx[ix][iy][iz]+=by[jx][jy][jz]*vx[kx][ky][kz]*DV[jx][jy][jz]; bvyy[ix][iy][iz]+=by[jx][jy][jz]*vy[kx][ky][kz]*DV[jx][jy][jz]; bvyz[ix][iy][iz]+=by[jx][jy][jz]*vz[kx][ky][kz]*DV[jx][jy][jz];
							bvzx[ix][iy][iz]+=bz[jx][jy][jz]*vx[kx][ky][kz]*DV[jx][jy][jz]; bvzy[ix][iy][iz]+=bz[jx][jy][jz]*vy[kx][ky][kz]*DV[jx][jy][jz]; bvzz[ix][iy][iz]+=bz[jx][jy][jz]*vz[kx][ky][kz]*DV[jx][jy][jz];
						 //vb
							vbxx[ix][iy][iz]+=vx[jx][jy][jz]*bx[kx][ky][kz]*DV[jx][jy][jz]; vbxy[ix][iy][iz]+=vx[jx][jy][jz]*by[kx][ky][kz]*DV[jx][jy][jz]; vbxz[ix][iy][iz]+=vx[jx][jy][jz]*bz[kx][ky][kz]*DV[jx][jy][jz];
							vbyx[ix][iy][iz]+=vy[jx][jy][jz]*bx[kx][ky][kz]*DV[jx][jy][jz]; vbyy[ix][iy][iz]+=vy[jx][jy][jz]*by[kx][ky][kz]*DV[jx][jy][jz]; vbyz[ix][iy][iz]+=vy[jx][jy][jz]*bz[kx][ky][kz]*DV[jx][jy][jz];
							vbzx[ix][iy][iz]+=vz[jx][jy][jz]*bx[kx][ky][kz]*DV[jx][jy][jz]; vbzy[ix][iy][iz]+=vz[jx][jy][jz]*by[kx][ky][kz]*DV[jx][jy][jz]; vbzz[ix][iy][iz]+=vz[jx][jy][jz]*bz[kx][ky][kz]*DV[jx][jy][jz];
						}
					}
				}
			}
		}
	}

	for(int ix=0;ix<=Lx; ix++)
	{double Qx=-Qx0+ix*DQx;
		for(int iy=0;iy<=Ly; iy++)
		{double Qy=-Qy0+iy*DQy;
			for(int iz=0;iz<=Lz; iz++)
			{double Qz=-Qz0+iz*DQz;
			// velociped
				rnvx[ix][iy][iz]=(vvxx[ix][iy][iz]+bbxx[ix][iy][iz])*Qx+(vvxy[ix][iy][iz]+bbxy[ix][iy][iz])*Qy+(vvxz[ix][iy][iz]+bbxz[ix][iy][iz])*Qz;
				rnvy[ix][iy][iz]=(vvyx[ix][iy][iz]+bbyx[ix][iy][iz])*Qx+(vvyy[ix][iy][iz]+bbyy[ix][iy][iz])*Qy+(vvyz[ix][iy][iz]+bbyz[ix][iy][iz])*Qz;
				rnvz[ix][iy][iz]=(vvzx[ix][iy][iz]+bbzx[ix][iy][iz])*Qx+(vvzy[ix][iy][iz]+bbzy[ix][iy][iz])*Qy+(vvzz[ix][iy][iz]+bbzz[ix][iy][iz])*Qz;
			// magnitno pole
				rnbx[ix][iy][iz]=(bvxx[ix][iy][iz]-vbxx[ix][iy][iz])*Qx+(bvxy[ix][iy][iz]-vbxy[ix][iy][iz])*Qy+(bvxz[ix][iy][iz]-vbxz[ix][iy][iz])*Qz;
				rnby[ix][iy][iz]=(bvyx[ix][iy][iz]-vbyx[ix][iy][iz])*Qx+(bvyy[ix][iy][iz]-vbyy[ix][iy][iz])*Qy+(bvyz[ix][iy][iz]-vbyz[ix][iy][iz])*Qz;
				rnbz[ix][iy][iz]=(bvzx[ix][iy][iz]-vbzx[ix][iy][iz])*Qx+(bvzy[ix][iy][iz]-vbzy[ix][iy][iz])*Qy+(bvzz[ix][iy][iz]-vbzz[ix][iy][iz])*Qz;
			}
		}
	}

	for(int ix=0;ix<=Lx; ix++)
	{double Qx=-Qx0+ix*DQx;
		for(int iy=0;iy<=Ly; iy++)
		{double Qy=-Qy0+iy*DQy;
			for(int iz=0;iz<=Lz; iz++)
			{double Qz=-Qz0+iz*DQz;
		     Q=sqrt(Qx*Qx+Qy*Qy+Qz*Qz); if (Q==0.0) continue;
			 rnx=Qx/Q; rny= Qy/Q; rnz= Qz/Q;
			 vn=rnvx[ix][iy][iz]*rnx+ rnvy[ix][iy][iz]*rny+rnvz[ix][iy][iz]*rnz;
			 bn=rnbx[ix][iy][iz]*rnx+ rnby[ix][iy][iz]*rny+rnbz[ix][iy][iz]*rnz; // bn ==? 0.
			 // velociped projection
			 rnvx[ix][iy][iz]-=rnx*vn; rnvy[ix][iy][iz]-=rny*vn; rnvz[ix][iy][iz]-=rnz*vn;
			 // magnetic projection
			 rnbx[ix][iy][iz]-=rnx*bn; rnby[ix][iy][iz]-=rny*bn; rnbz[ix][iy][iz]-=rnz*bn;
			 printf("%g %g %g %g \n",Qx, Qy, Qz, vvxx[ix][iy][iz] );
			}
		}
	}
	// vvxx[ix][iy][iz] - exp(-.5*(Qx*Qx+Qy*Qy+Qz*Qz))/(sqrt(8.*M_PI)*sqrt(8.*M_PI)*sqrt(8.*M_PI));
	
	
	double alpha=0.0;
double betax,betaz;
double betay=0.0; 
double omega=-2./3.;
double fnvx,fnvy,fnvz,fnbx,fnby,fnbz;

betax = sin(alpha); 
betaz = cos(alpha);


	for(int ix=0;ix<=Lx; ix++)
	{double Qx=-Qx0+ix*DQx;
		for(int iy=0;iy<=Ly; iy++)
		{double Qy=-Qy0+iy*DQy;
			for(int iz=0;iz<=Lz; iz++)
			{double Qz=-Qz0+iz*DQz;

             Q2=Qx*Qx+Qy*Qy+Qz*Qz;
			 Q=sqrt(Q2)
			 rnx=Qx/Q; rny=Qy/Q; rnz=Qz/Q;
			 Qbeta=Qx*betax+Qz*betaz; // betay=0.
			 					 
			 fnvx = (vvxx[ix][iy][iz] + bbxx[ix][iy][iz] )* Qx + (vvxy[ix][iy][iz] + bbxy[ix][iy][iz]) * Qy + (vvxz[ix][iy][iz] + bbxz[ix][iy][iz]) * Qz;
			 fnvy = (vvyx[ix][iy][iz] + bbyx[ix][iy][iz] )* Qx + (vvyy[ix][iy][iz] + bbyy[ix][iy][iz]) * Qy + (vvyz[ix][iy][iz] + bbyz[ix][iy][iz]) * Qz;
			 fnvz = (vvzx[ix][iy][iz] + bbzx[ix][iy][iz] )* Qx + (vvzy[ix][iy][iz] + bbzy[ix][iy][iz]) * Qy + (vvzz[ix][iy][iz] + bbzz[ix][iy][iz]) * Qz;

			 fnbx = (bvxx[ix][iy][iz] - vbxx[ix][iy][iz] )* Qx + (bvxy[ix][iy][iz] - vbxy[ix][iy][iz]) * Qy + (bvxz[ix][iy][iz] - vbxz[ix][iy][iz]) * Qz;
			 fnby = (bvyx[ix][iy][iz] - vbyx[ix][iy][iz] )* Qx + (bvyy[ix][iy][iz] - vbyy[ix][iy][iz]) * Qy + (bvyz[ix][iy][iz] - vbyz[ix][iy][iz]) * Qz;
			 fnbz = (bvzx[ix][iy][iz] - vbzx[ix][iy][iz] )* Qx + (bvzy[ix][iy][iz] - vbzy[ix][iy][iz]) * Qy + (bvzz[ix][iy][iz] - vbzz[ix][iy][iz]) * Qz;

			 fvn = fnvx*rnx + fnvy*rny + fnvz*rnz; 
			 fbn = fnbx*rnx + fnby*rny + fnbz*rnz;

			 fnvx-=fvn*rnx;
			 fnvy-=fvn*rny;
			 fnvz-=fvn*rnz;

			 fnbx-=fbn*rnx;
			 fnby-=fbn*rny;
			 fnbz-=fbn*rnz;

			 aaa=2.0*(ny*vx[ix][iy][iz] - omega*(nx*vy[ix][iy][iz]-ny*vx[ix][iy][iz]));
			 flvx=                  aaa*rnx - 2.0*omega*vy[ix][iy][iz] + Qbeta*bx[ix][iy][iz] - rnuk*Q2*vx[ix][iy][iz];
			 flvy=-vx[ix][iy][iz] + aaa*rny + 2.0*omega*vx[ix][iy][iz] + Qbeta*by[ix][iy][iz] - rnuk*Q2*vy[ix][iy][iz];
			 flvz=                  aaa*rnz                            + Qbeta*bz[ix][iy][iz] - rnuk*Q2*vz[ix][iy][iz];

			 flbx=      aaa*rnx - 2.0*omega*vy[ix][iy][iz] - Qbeta*vx[ix][iy][iz] - rnum*Q2*vx[ix][iy][iz];
			 flby=+bx + aaa*rny + 2.0*omega*vx[ix][iy][iz] - Qbeta*vy[ix][iy][iz] - rnum*Q2*vy[ix][iy][iz];
			 flbz=      aaa*rnz                            - Qbeta*vz[ix][iy][iz] - rnum*Q2*vz[ix][iy][iz];

			 fvx[ix][iy][iz]=flvx+fnvx; 
			 fvy[ix][iy][iz]=flvy+fnvy; 
			 fvz[ix][iy][iz]=flvz+fnvz; 

			 fbx[ix][iy][iz]=flbx+fnbx; 
			 fby[ix][iy][iz]=flby+fnby; 
			 fbz[ix][iy][iz]=flbz+fnbz;
			 
			}
		}
	}

	return 0;
}

