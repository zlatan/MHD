#include "parameters.h"
#include "variables.h"
#include <stdlib.h>
#include <math.h>


void calculate_nonlinear_term(NonlinearTerm *bb, NonlinearTerm *vv, NonlinearTerm *bv, NonlinearTerm *vb,LinearTerm *b, LinearTerm *v, double DQx,double DQy,double DQz,double DVQ, double DV[2*Nx][2*Nx][2*Nx])
{
		   for(int ix=0;ix<=2*Nx; ix++)
		       {double Qx=-Qx0+ix*DQx;
		       for(int iy=0;iy<=2*Ny; iy++)
			   {double Qy=-Qy0+iy*DQy;
			   for(int iz=0;iz<=2*Nz; iz++)
			       {double Qz=-Qz0+iz*DQz;
			       // summation cycles for nonlinear terms
			       for (int jx=0;jx<=2*Nx; jx++)
				   {int kx=ix-jx+Nx; 
				    if(kx<0) continue;
                                    if(kx>2*Nx) continue;
				    for (int jy=0;jy<=2*Ny; jy++)
					{int ky=iy-jy+Ny; 
                                         if(ky<0) continue;
                                         if(ky>2*Ny) continue;
					 for (int jz=0;jz<=2*Nz; jz++)
					 {int kz=iz-jz+Nz; if(kz<0) continue; if(kz>2*Nz) continue;					   
						//bb
						bb->xx[ix][iy][iz]+=b->x[jx][jy][jz]*b->x[kx][ky][kz]*DV[jx][jy][jz];
						bb->xy[ix][iy][iz]+=b->x[jx][jy][jz]*b->y[kx][ky][kz]*DV[jx][jy][jz];
 						bb->xz[ix][iy][iz]+=b->x[jx][jy][jz]*b->z[kx][ky][kz]*DV[jx][jy][jz];

						bb->yx[ix][iy][iz]+=b->y[jx][jy][jz]*b->x[kx][ky][kz]*DV[jx][jy][jz];
						bb->yy[ix][iy][iz]+=b->y[jx][jy][jz]*b->y[kx][ky][kz]*DV[jx][jy][jz];
						bb->yz[ix][iy][iz]+=b->y[jx][jy][jz]*b->z[kx][ky][kz]*DV[jx][jy][jz];

						bb->zx[ix][iy][iz]+=b->z[jx][jy][jz]*b->x[kx][ky][kz]*DV[jx][jy][jz];
						bb->zy[ix][iy][iz]+=b->z[jx][jy][jz]*b->y[kx][ky][kz]*DV[jx][jy][jz];
 						bb->zz[ix][iy][iz]+=b->z[jx][jy][jz]*b->z[kx][ky][kz]*DV[jx][jy][jz];
						//vv
						vv->xx[ix][iy][iz]+=v->x[jx][jy][jz]*v->x[kx][ky][kz]*DV[jx][jy][jz];
						vv->xy[ix][iy][iz]+=v->x[jx][jy][jz]*v->y[kx][ky][kz]*DV[jx][jy][jz];
		 				vv->xz[ix][iy][iz]+=v->x[jx][jy][jz]*v->z[kx][ky][kz]*DV[jx][jy][jz];

						vv->yx[ix][iy][iz]+=v->y[jx][jy][jz]*v->x[kx][ky][kz]*DV[jx][jy][jz];
 						vv->yy[ix][iy][iz]+=v->y[jx][jy][jz]*v->y[kx][ky][kz]*DV[jx][jy][jz]; 
						vv->yz[ix][iy][iz]+=v->y[jx][jy][jz]*v->z[kx][ky][kz]*DV[jx][jy][jz];
						
						vv->zx[ix][iy][iz]+=v->z[jx][jy][jz]*v->x[kx][ky][kz]*DV[jx][jy][jz];
 						vv->zy[ix][iy][iz]+=v->z[jx][jy][jz]*v->y[kx][ky][kz]*DV[jx][jy][jz];
						vv->zz[ix][iy][iz]+=v->z[jx][jy][jz]*v->z[kx][ky][kz]*DV[jx][jy][jz];
						//bv
						bv->xx[ix][iy][iz]+=b->x[jx][jy][jz]*v->x[kx][ky][kz]*DV[jx][jy][jz];
						bv->xy[ix][iy][iz]+=b->x[jx][jy][jz]*v->y[kx][ky][kz]*DV[jx][jy][jz];
						bv->xz[ix][iy][iz]+=b->x[jx][jy][jz]*v->z[kx][ky][kz]*DV[jx][jy][jz];

						bv->yx[ix][iy][iz]+=b->y[jx][jy][jz]*v->x[kx][ky][kz]*DV[jx][jy][jz];
						bv->yy[ix][iy][iz]+=b->y[jx][jy][jz]*v->y[kx][ky][kz]*DV[jx][jy][jz];
						bv->yz[ix][iy][iz]+=b->y[jx][jy][jz]*v->z[kx][ky][kz]*DV[jx][jy][jz];

						bv->zx[ix][iy][iz]+=b->z[jx][jy][jz]*v->x[kx][ky][kz]*DV[jx][jy][jz];
						bv->zy[ix][iy][iz]+=b->z[jx][jy][jz]*v->y[kx][ky][kz]*DV[jx][jy][jz];
						bv->zz[ix][iy][iz]+=b->z[jx][jy][jz]*v->z[kx][ky][kz]*DV[jx][jy][jz];
						//vb
						vb->xx[ix][iy][iz]+=v->x[jx][jy][jz]*b->x[kx][ky][kz]*DV[jx][jy][jz];
						vb->xy[ix][iy][iz]+=v->x[jx][jy][jz]*b->y[kx][ky][kz]*DV[jx][jy][jz];
						vb->xz[ix][iy][iz]+=v->x[jx][jy][jz]*b->z[kx][ky][kz]*DV[jx][jy][jz];
												
						vb->yx[ix][iy][iz]+=v->y[jx][jy][jz]*b->x[kx][ky][kz]*DV[jx][jy][jz];
 						vb->yy[ix][iy][iz]+=v->y[jx][jy][jz]*b->y[kx][ky][kz]*DV[jx][jy][jz];
						vb->yz[ix][iy][iz]+=v->y[jx][jy][jz]*b->z[kx][ky][kz]*DV[jx][jy][jz];

						vb->zx[ix][iy][iz]+=v->z[jx][jy][jz]*b->x[kx][ky][kz]*DV[jx][jy][jz];
						vb->zy[ix][iy][iz]+=v->z[jx][jy][jz]*b->y[kx][ky][kz]*DV[jx][jy][jz];
						vb->zz[ix][iy][iz]+=v->z[jx][jy][jz]*b->z[kx][ky][kz]*DV[jx][jy][jz];
						}//kz
					    }//ky
				        }//kx
			           }//iz
			      }//iy
		         }//ix
}