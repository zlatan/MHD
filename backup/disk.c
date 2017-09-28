#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void)
{

	double Qx0=20., Qy0=3., Qz0=3.;
	Qx0=Qy0=Qz0=.;// xxxx
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
	double rnuk=0.01;
	double rnum=0.01;
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
		        DV[ix][iy][iz]=DVQ;
/*
				bbxx[ix][iy][iz]=(rand()-.5)*1E-10;
				bbxy[ix][iy][iz]=(rand()-.5)*1E-10;
				bbxz[ix][iy][iz]=(rand()-.5)*1E-10;
				bbyx[ix][iy][iz]=(rand()-.5)*1E-10;
				bbyy[ix][iy][iz]=(rand()-.5)*1E-10;
				bbyz[ix][iy][iz]=(rand()-.5)*1E-10;
				bbzx[ix][iy][iz]=(rand()-.5)*1E-10;
				bbzy[ix][iy][iz]=(rand()-.5)*1E-10;
				bbzz[ix][iy][iz]=(rand()-.5)*1E-10;
				vvxx[ix][iy][iz]=(rand()-.5)*1E-10;
				vvxy[ix][iy][iz]=(rand()-.5)*1E-10;
				vvxz[ix][iy][iz]=(rand()-.5)*1E-10;
				vvyx[ix][iy][iz]=(rand()-.5)*1E-10;
				vvyy[ix][iy][iz]=(rand()-.5)*1E-10;
				vvyz[ix][iy][iz]=(rand()-.5)*1E-10;
				vvzx[ix][iy][iz]=(rand()-.5)*1E-10;
				vvzy[ix][iy][iz]=(rand()-.5)*1E-10;
				vvzz[ix][iy][iz]=(rand()-.5)*1E-10;
				bvxx[ix][iy][iz]=(rand()-.5)*1E-10;
				bvxy[ix][iy][iz]=(rand()-.5)*1E-10;
				bvxz[ix][iy][iz]=(rand()-.5)*1E-10;
				bvyx[ix][iy][iz]=(rand()-.5)*1E-10;
				bvyy[ix][iy][iz]=(rand()-.5)*1E-10;
				bvyz[ix][iy][iz]=(rand()-.5)*1E-10;
				bvzx[ix][iy][iz]=(rand()-.5)*1E-10;
				bvzy[ix][iy][iz]=(rand()-.5)*1E-10;
				bvzz[ix][iy][iz]=(rand()-.5)*1E-10;
				vbxx[ix][iy][iz]=(rand()-.5)*1E-10;
				vbxy[ix][iy][iz]=(rand()-.5)*1E-10;
				vbxz[ix][iy][iz]=(rand()-.5)*1E-10;
				vbyx[ix][iy][iz]=(rand()-.5)*1E-10;
				vbyy[ix][iy][iz]=(rand()-.5)*1E-10;
				vbyz[ix][iy][iz]=(rand()-.5)*1E-10;
				vbzx[ix][iy][iz]=(rand()-.5)*1E-10;
				vbzy[ix][iy][iz]=(rand()-.5)*1E-10;
				vbzz[ix][iy][iz]=(rand()-.5)*1E-10;
*/

				vx[ix][iy][iz] = 1.0;
				vy[ix][iy][iz] = 1.0;
				vz[ix][iy][iz] = 0.0;

				bx[ix][iy][iz] = 1.0;
				by[ix][iy][iz] = 1.0;
				bz[ix][iy][iz] = 0.0;


				// vx[ix][iy][iz]=exp(-(Qx*Qx+Qy*Qy+Qz*Qz));// inicializacija na skorosti

				Q2=Qx*Qx+Qy*Qy+Qz*Qz;
				Q=sqrt(Q2);
				// printf("Q=%g\n",Q);
				if ( ((ix - Nx)*(ix - Nx) + (iy-Ny)*(iy-Ny) + (iz-Nz)*(iz-Nz)) == 0 )
				{
					// printf("Zero -> Q=%g\n",Q);
					// printf("ix-Nx=%i iy-Ny=%i iz-Nz=%i\n",ix-Nx,iy-Ny,iz-Nz);
					continue;
				}
				rnx=Qx/Q; rny=Qy/Q; rnz=Qz/Q;

				vn = vx[ix][iy][iz]*rnx + vy[ix][iy][iz]*rny + vz[ix][iy][iz]*rnz;
				bn = bx[ix][iy][iz]*rnx + by[ix][iy][iz]*rny + bz[ix][iy][iz]*rnz;

				vx[ix][iy][iz]-= vn*rnx;
				vy[ix][iy][iz]-= vn*rny;
				vz[ix][iy][iz]-= vn*rnz;

				bx[ix][iy][iz]-= bn*rnx;
				by[ix][iy][iz]-= bn*rny;
				bz[ix][iy][iz]-= bn*rnz;
				// printf("%g %g %g %g %g %g\n",vvxx[ix][iy][iz],bbxx[ix][iy][iz],vvxy[ix][iy][iz],bbxy[ix][iy][iz],vvxz[ix][iy][iz],bbxz[ix][iy][iz]);

			}
		}
	}

		Tempi=10.; Nt=2; dt=Tempi/Nt;
        // nachalo na integrirane po vreme; tuk se vyrti
		for(it=1;it<=Nt; it++)
		{t=it*dt;
		  for (int ix=0;ix<=Lx; ix++)
				{double Qx=-Qx0+ix*DQx;
				for (int iy=0;iy<=Ly; iy++)
				{double Qy=-Qy0+iy*DQy;
						for (int iz=0;iz<=Lz; iz++)
								{double Qz=-Qz0+iz*DQz;
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

												// printf("%g %g %g %g %g %g\n",vvxx[ix][iy][iz],bbxx[ix][iy][iz],vvxy[ix][iy][iz],bbxy[ix][iy][iz],vvxz[ix][iy][iz],bbxz[ix][iy][iz]);

											}//kz
										}//ky
									}//kx
								}//iz
							}//iy
						}//ix

		 // krai na nelinejnata sila

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
						 printf("Zero");
						 continue;
					 }
					 Q=sqrt(Q2);
					//  printf("Q=%g it=>%i",Q,it);
					 rnx=Qx/Q; rny=Qy/Q; rnz=Qz/Q;
					 Qbeta=Qx*betax+Qz*betaz; // betay=0.

					 fnvx = (vvxx[ix][iy][iz] + bbxx[ix][iy][iz] )* Qx + (vvxy[ix][iy][iz] + bbxy[ix][iy][iz]) * Qy + (vvxz[ix][iy][iz] + bbxz[ix][iy][iz]) * Qz;
					 fnvy = (vvyx[ix][iy][iz] + bbyx[ix][iy][iz] )* Qx + (vvyy[ix][iy][iz] + bbyy[ix][iy][iz]) * Qy + (vvyz[ix][iy][iz] + bbyz[ix][iy][iz]) * Qz;
					 fnvz = (vvzx[ix][iy][iz] + bbzx[ix][iy][iz] )* Qx + (vvzy[ix][iy][iz] + bbzy[ix][iy][iz]) * Qy + (vvzz[ix][iy][iz] + bbzz[ix][iy][iz]) * Qz;

					//  printf("fnvx=%g fnvy=%g fnvz=%g iteration=>%i\n",fnvx,fnvy,fnvz,it);
					//  printf("rnx=%g rny=%g rnz=%g iteration=>%i\n",rnx,rny,rnz,it);

					 fnbx = (bvxx[ix][iy][iz] - vbxx[ix][iy][iz] )* Qx + (bvxy[ix][iy][iz] - vbxy[ix][iy][iz]) * Qy + (bvxz[ix][iy][iz] - vbxz[ix][iy][iz]) * Qz;
					 fnby = (bvyx[ix][iy][iz] - vbyx[ix][iy][iz] )* Qx + (bvyy[ix][iy][iz] - vbyy[ix][iy][iz]) * Qy + (bvyz[ix][iy][iz] - vbyz[ix][iy][iz]) * Qz;
					 fnbz = (bvzx[ix][iy][iz] - vbzx[ix][iy][iz] )* Qx + (bvzy[ix][iy][iz] - vbzy[ix][iy][iz]) * Qy + (bvzz[ix][iy][iz] - vbzz[ix][iy][iz]) * Qz;

					 fvn = fnvx*rnx + fnvy*rny + fnvz*rnz;
					 fbn = fnbx*rnx + fnby*rny + fnbz*rnz;

					//  printf("fvn=%g it=>%i",fvn,it);

					 fnvx-= fvn*rnx;
					 fnvy-= fvn*rny;
					 fnvz-= fvn*rnz;

					//  printf("fnvx=%g fnvy=%g fnvz=%g iteration=>%i\n",fnvx,fnvy,fnvz,it);

					 fnbx-= fbn*rnx;
					 fnby-= fbn*rny;
					 fnbz-= fbn*rnz;

					 aaa = 2.0*(rny*vx[ix][iy][iz] - omega*(rnx*vy[ix][iy][iz]-rny*vx[ix][iy][iz]));
					 flvx =                  aaa*rnx + 2.0*omega*vy[ix][iy][iz] + Qbeta*bx[ix][iy][iz] - rnuk*Q2*vx[ix][iy][iz];
					 flvy =-vx[ix][iy][iz] + aaa*rny - 2.0*omega*vx[ix][iy][iz] + Qbeta*by[ix][iy][iz] - rnuk*Q2*vy[ix][iy][iz];
					 flvz =                  aaa*rnz                            + Qbeta*bz[ix][iy][iz] - rnuk*Q2*vz[ix][iy][iz];

					 flbx =    		- Qbeta*vx[ix][iy][iz] - rnum*Q2*vx[ix][iy][iz];
					 flby = bx[ix][iy][iz]  - Qbeta*vy[ix][iy][iz] - rnum*Q2*vy[ix][iy][iz];
					 flbz =    		- Qbeta*vz[ix][iy][iz] - rnum*Q2*vz[ix][iy][iz];


					 fvn = flvx*rnx + flvy*rny + flvz*rnz;
					 fbn = flbx*rnx + flby*rny + flbz*rnz;

					//  printf("fvn=%g fbn=%g iteration=>%i\n",fvn,fbn,it);

					 flvx-= fvn*rnx;
					 flvy-= fvn*rny;
					 flvz-= fvn*rnz;

					 flbx-= fbn*rnx;
					 flby-= fbn*rny;
					 flbz-= fbn*rnz;

					//  printf("flvx=%g flvy=%g flvz=%g iteration=>%i\n",flvx,flvy,flvz,it);


					 fvx[ix][iy][iz] = flvx+fnvx;
					 fvy[ix][iy][iz] = flvy+fnvy;
					 fvz[ix][iy][iz] = flvz+fnvz;

					//  printf("x->%g y->%g z->%g \n",fvx[ix][iy][iz],fvy[ix][iy][iz],fvy[ix][iy][iz]);

					 fbx[ix][iy][iz] = flbx+fnbx;
					 fby[ix][iy][iz] = flby+fnby;
					 fbz[ix][iy][iz] = flbz+fnbz;

					 vx[ix][iy][iz]+=fvx[ix][iy][iz]*dt; vy[ix][iy][iz]+=fvy[ix][iy][iz]*dt; vz[ix][iy][iz]+=fvz[ix][iy][iz]*dt;
           bx[ix][iy][iz]+=fbx[ix][iy][iz]*dt; by[ix][iy][iz]+=fby[ix][iy][iz]*dt; bz[ix][iy][iz]+=fbz[ix][iy][iz]*dt;

					//  printf("f-> %g %g %g %i\n",fbx[ix][iy][iz], fby[ix][iy][iz], fbz[ix][iy][iz],it);
						if (it == 1)
						{
							Ei[ix][iy][iz] = vx[ix][iy][iz]*vx[ix][iy][iz] + vy[ix][iy][iz]*vy[ix][iy][iz] + vz[ix][iy][iz]*vz[ix][iy][iz] + bx[ix][iy][iz]*bx[ix][iy][iz] + by[ix][iy][iz]*by[ix][iy][iz] + bz[ix][iy][iz]*bz[ix][iy][iz];
						}

						if (it == 2)
						{
							Ef[ix][iy][iz] = vx[ix][iy][iz]*vx[ix][iy][iz] + vy[ix][iy][iz]*vy[ix][iy][iz] + vz[ix][iy][iz]*vz[ix][iy][iz] + bx[ix][iy][iz]*bx[ix][iy][iz] + by[ix][iy][iz]*by[ix][iy][iz] + bz[ix][iy][iz]*bz[ix][iy][iz];
						}
					//  printf("(v^2(Q) + b^2(Q))[t]-> %g %g %i\n",vx[ix][iy][iz]*vx[ix][iy][iz] + vy[ix][iy][iz]*vy[ix][iy][iz] + vz[ix][iy][iz]*vz[ix][iy][iz] + bx[ix][iy][iz]*bx[ix][iy][iz] + by[ix][iy][iz]*by[ix][iy][iz] + bz[ix][iy][iz]*bz[ix][iy][iz],Q,it);
					}//iz
				}//iy
			}//ix
		}//it


					for(int ix=0;ix<=Lx; ix++)
					{double Qx=-Qx0+ix*DQx;
					// 	for(int iy=0;iy<=Ly; iy++)
					// 	{double Qy=-Qy0+iy*DQy;
							for(int iz=0;iz<=Lz; iz++)
							{double Qz=-Qz0+iz*DQz;
								// if ( ((ix - Nx)*(ix - Nx) + (iy-Ny)*(iy-Ny) + (iz-Nz)*(iz-Nz)) == 0 )
								// {
								// 	// printf("Zero");
								// 	continue;
								// }
								printf("%g %g %g\n",Qx,Qz,pow( (1/(2.0*Tempi))*log(Ef[ix][Ny][iz]/Ei[ix][Ny][iz]),2));
 
							}//iz
					// 	}//iy
					}//ix





	return 0;
}// main
