#incude <stdio.h>

int main(void)
{
	
	double Qx0=20., Qy0=Qz0=3., DQx, DQy, Dqz;
	int Nx=30, Ny=Nz=6;
	int NDimx=2*Nx+2, NDimy=2*Ny+2, NDimz=2*Nz+2;
	double vx[NDimx][NDimy][NDimz], vy[NDimx][NDimy][NDimz], vz[NDimx][NDimy][NDimz], 
	bx[NDimx][NDimy][NDimz], by[NDimx][NDimy][NDimz], bz[NDimx][NDimy][NDimz], 
	bbxx[NDimx][NDimy][NDimz], bbxy[NDimx][NDimy][NDimz], bbxz[NDimx][NDimy][NDimz], 
	bbyx[NDimx][NDimy][NDimz], bbyy[NDimx][NDimy][NDimz], bbyz[NDimx][NDimy][NDimz], 
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
	vx[NDimx][NDimy][NDimz], vy[NDimx][NDimy][NDimz], vz[NDimx][NDimy][NDimz], 
	rnvx[NDimx][NDimy][NDimz], rnvy[NDimx][NDimy][NDimz], rnvz[NDimx][NDimy][NDimz], 
	rnbx[NDimx][NDimy][NDimz], rnby[NDimx][NDimy][NDimz], rnbz[NDimx][NDimy][NDimz]; 
	int Lx=2*Nx, Ly= 2*Ny, Lz=2*Nz;
	double DQx=2*Qx0/(2*Nx), DQy=2*Qy0/(2*Ny), DQx=2*Qz0/(2*Nz), rnx, rny, rnz, vn, bn;
	
	for (ix=0;ix<=Lx; ix++)
	{//Qx=-Qx0+ix*DQx;
		for (iy=0;iy<=Ly; iy++)	
		{//Qy=-Qy0+iy*DQy;
			for (iz=0;iz<=Lz; iz++)
			{//Qz=-Qz0+iz*DQz;
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
		vbzx[ix][iy][iz]= vbzy[ix][iy][iz]= vbzz[ix][iy][iz]=
	    0.0;
			// summation cycles for nonlinear terms
				for (jx=0;jx<=Lx; jx++)
				{kx=ix-jx+Nx; if(kx<0) continue; if(kx>Lx) continue; 
					for (jy=0;jy<=Ly; jy++)
					{ky=iy-jy+Ny; if(ky<0) continue; if(ky>Ly) continue; 
						for (jz=0;jz<=Lz; jz++)
						{// razgele, naj setne zapoichvame rabota
						 kz=iz-jz+Nz; if(kz<0) continue; if(kz>Lz) continue;
						 //bb
					     	bbxx[ix][iy][iz]+=bx[jx][jy][jz]*bx[kx][ky][kz]; bbxy[ix][iy][iz]+=bx[jx][jy][jz]*by[kx][ky][kz]; bbxz[ix][iy][iz]+=bx[jx][jy][jz]*bz[kx][ky][kz];
							bbyx[ix][iy][iz]+=by[jx][jy][jz]*bx[kx][ky][kz]; bbyy[ix][iy][iz]+=by[jx][jy][jz]*by[kx][ky][kz]; bbyz[ix][iy][iz]+=by[jx][jy][jz]*bz[kx][ky][kz];
							bbzx[ix][iy][iz]+=bz[jx][jy][jz]*bx[kx][ky][kz]; bbzy[ix][iy][iz]+=by[jx][jy][jz]*by[kx][ky][kz]; bbzz[ix][iy][iz]+=bz[jx][jy][jz]*bz[kx][ky][kz];
						 //vv	
							vvxx[ix][iy][iz]+=vx[jx][jy][jz]*vx[kx][ky][kz]; vvxy[ix][iy][iz]+=vx[jx][jy][jz]*vy[kx][ky][kz]; vvxz[ix][iy][iz]+=vx[jx][jy][jz]*vz[kx][ky][kz];
							vvyx[ix][iy][iz]+=vy[jx][jy][jz]*vx[kx][ky][kz]; vvyy[ix][iy][iz]+=vy[jx][jy][jz]*vy[kx][ky][kz]; vvyz[ix][iy][iz]+=vy[jx][jy][jz]*vz[kx][ky][kz];
							vvzx[ix][iy][iz]+=vz[jx][jy][jz]*vx[kx][ky][kz]; vvzy[ix][iy][iz]+=vy[jx][jy][jz]*vy[kx][ky][kz]; vvzz[ix][iy][iz]+=vz[jx][jy][jz]*vz[kx][ky][kz];
						 //bv	
							bvxx[ix][iy][iz]+=bx[jx][jy][jz]*vx[kx][ky][kz]; bvxy[ix][iy][iz]+=bx[jx][jy][jz]*vy[kx][ky][kz]; bvxz[ix][iy][iz]+=bx[jx][jy][jz]*vv[kx][ky][kz];
							bvyx[ix][iy][iz]+=by[jx][jy][jz]*vx[kx][ky][kz]; bbyy[ix][iy][iz]+=by[jx][jy][jz]*vy[kx][ky][kz]; bbyz[ix][iy][iz]+=by[jx][jy][jz]*vz[kx][ky][kz];
							bvzx[ix][iy][iz]+=bz[jx][jy][jz]*vx[kx][ky][kz]; bbzy[ix][iy][iz]+=by[jx][jy][jz]*vy[kx][ky][kz]; bbzz[ix][iy][iz]+=bz[jx][jy][jz]*vz[kx][ky][kz];
						 //vb
						    vvxx[ix][iy][iz]+=vx[jx][jy][jz]*bx[kx][ky][kz]; vvxy[ix][iy][iz]+=vx[jx][jy][jz]*by[kx][ky][kz]; vvxz[ix][iy][iz]+=vx[jx][jy][jz]*bz[kx][ky][kz];
							vvyx[ix][iy][iz]+=vy[jx][jy][jz]*bx[kx][ky][kz]; vvyy[ix][iy][iz]+=vy[jx][jy][jz]*by[kx][ky][kz]; vvyz[ix][iy][iz]+=vy[jx][jy][jz]*bz[kx][ky][kz];
							vvzx[ix][iy][iz]+=vz[jx][jy][jz]*bx[kx][ky][kz]; vvzy[ix][iy][iz]+=vy[jx][jy][jz]*by[kx][ky][kz]; vvzz[ix][iy][iz]+=vz[jx][jy][jz]*bz[kx][ky][kz];
						}	
					}
				}					
			}	
		}	
	}	
	
	for (ix=0;ix<=Lx; ix++)
	{Qx=-Qx0+ix*DQx;
		for (iy=0;iy<=Ly; iy++)	
		{Qy=-Qy0+iy*DQy;
			for (iz=0;iz<=Lz; iz++)
			{Qz=-Qz0+iz*DQz;
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
	
	for (ix=0;ix<=Lx; ix++)
	{Qx=-Qx0+ix*DQx;
		for (iy=0;iy<=Ly; iy++)	
		{Qy=-Qy0+iy*DQy;	
			for (iz=0;iz<=Lz; iz++)
			{Qz=-Qz0+iz*DQz;
		     Q=sqrt(Qx*Qx+Qy*Qy+Qz*Qz); if (Q==0.0) continue;
			 rnx=Qx/Q; rny= Qy/Q; rnz= Qz/Q;
			 vn=rnvx[ix][iy][iz]*rnx+ rnvy[ix][iy][iz]*rny+rnvz[ix][iy][iz]*rnz;
			 bn=rnbx[ix][iy][iz]*rnx+ rnby[ix][iy][iz]*rny+rnbz[ix][iy][iz]*rnz; // bn ==? 0.
			 // velociped projection 
			 rnvx[ix][iy][iz]-=rnx*vn; rnvy[ix][iy][iz]-=rny*vn; rnvz[ix][iy][iz]-=rnz*vn;
			 // magnetic projection 
			 rnbx[ix][iy][iz]-=rnb*bn; rnby[ix][iy][iz]-=rny*bn; rnbz[ix][iy][iz]-=rnz*bn;
			}	
		}
	}
	
	return 0;
}
