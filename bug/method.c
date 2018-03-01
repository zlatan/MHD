#include <stdio.h>
#include <math.h>

double w=-2.0/3.0;

void point(double Qx, double Qy, double Qz, double vx, double vy, double vz, double bx, double by, double bz)
{

double Q[3];
double v[3],b[3];
double rn[3];
double fv[3];
double fb[3];

double vv[3][3];
double vb[3][3];
double bv[3][3];
double bb[3][3];
double fvl=0.0;
double fbl=0.0;


                                         for(int i=1;i<=3;i++) 
                                           for(int j=1;j<=3;j++)  
                                           {
                                             vv[i][j]=0;
                                             vb[i][j]=0;
                                             bv[i][j]=0;
                                             bb[i][j]=0;
                                           }  

double Q2,Qalpha;



                                         // begin Todor Method Area
                                         // -----------------------------------------------------------------------------   
                                         // nachalo na prisvojavane na lokalni promenlivi
                                         Q[1]=Qx;  Q[2]=Qy;  Q[3]=Qz; 
/*
                                         v[1]=vx[ix][iy][iz]; v[2]=vy[ix][iy][iz]; v[3]=vz[ix][iy][iz];
                                         b[1]=bx[ix][iy][iz]; b[2]=vy[ix][iy][iz]; b[3]=vz[ix][iy][iz];
*/

                                         v[1]=vx; v[2]=vy; v[3]=vz;
                                         b[1]=bx; b[2]=vy; b[3]=vz;

/*
	                                 vv[1][1]=wwxx[ix][iy][iz]; vv[1][2]=wwxy[ix][iy][iz]; vv[1][3]=wwxz[ix][iy][iz];
                                         vv[2][1]=wwyx[ix][iy][iz]; vv[2][2]=wwyy[ix][iy][iz]; vv[2][3]=wwyz[ix][iy][iz];
                                         vv[3][1]=wwzx[ix][iy][iz]; vv[3][2]=wwzy[ix][iy][iz]; vv[3][3]=wwzz[ix][iy][iz];
                                         
                                         vb[1][1]=wwxx[ix][iy][iz]; vb[1][2]=wwxy[ix][iy][iz]; vb[1][3]=wwxz[ix][iy][iz];
                                         vb[2][1]=wwyx[ix][iy][iz]; vb[2][2]=wwyy[ix][iy][iz]; vb[2][3]=wwyz[ix][iy][iz];
                                         vb[3][1]=wwzx[ix][iy][iz]; vb[3][2]=wwzy[ix][iy][iz]; vb[3][3]=wwzz[ix][iy][iz];

                                         bv[1][1]=wwxx[ix][iy][iz]; bv[1][2]=wwxy[ix][iy][iz]; bv[1][3]=wwxz[ix][iy][iz];
                                         bv[2][1]=wwyx[ix][iy][iz]; bv[2][2]=wwyy[ix][iy][iz]; bv[2][3]=wwyz[ix][iy][iz];
                                         bv[3][1]=wwzx[ix][iy][iz]; bv[3][2]=wwzy[ix][iy][iz]; bv[3][3]=wwzz[ix][iy][iz];

                                         bb[1][1]=wwxx[ix][iy][iz]; bb[1][2]=wwxy[ix][iy][iz]; bb[1][3]=wwxz[ix][iy][iz];
                                         bb[2][1]=wwyx[ix][iy][iz]; bb[2][2]=wwyy[ix][iy][iz]; bb[2][3]=wwyz[ix][iy][iz];
                                         bb[3][1]=wwzx[ix][iy][iz]; bb[3][2]=wwzy[ix][iy][iz]; bb[3][3]=wwzz[ix][iy][iz];
*/
                                         // kraj na prisvojavaneto na lokalni promenlivi
					 for(int i=1;i<=3;i++) 
					 {fv[i]=0.0;
					  fb[i]=0.0;} // nulirane na proizvodni
                                         fv[1]= 2.*w*v[2];
                                         fv[2]=-2.*w*v[1]; // Coriolis
					 for(int i=1;i<=3;i++) 
					 {fv[i]+= Qalpha*b[i]-rn[i]*Q2*v[i];
					  fb[i]+=-Qalpha*v[i]-rn[i]*Q2*b[i];} // ba li gy kvo znachi

                                         for(int i=1;i<=3;i++) 
					 {
                                           for(int j=1;j<=3;j++)  
                                           {
                                             fv[i]+=(vv[i][j]+vv[i][j])*Q[j];
                                             fb[i]+=(bv[i][j]-vb[i][j])*Q[j];
                                           }  
					 } // nelinejni chlenove


					 for(int i=1;i<=3;i++) 
					 {fvl+=rn[i]*fv[i];
					  fbl+=rn[i]*fb[i];} // nadlyzhni chlenove
                                         for(int i=1;i<=3;i++) 
					 {fv[i]-=rn[i]*fvl;
					  fb[i]-=rn[i]*fbl;} // ortogonalizacija
                                         
                                         for(int i=1;i<=3;i++) 
					 {fv[i]+=2.0*rn[2]*v[1]*rn[i];} // srjazvane
					 fv[2]-=v[1];
                                         fb[2]+=b[1]; // kraj na srqzvaneto; ednotipni clenove

					
                                         for(int i=1;i<=3;i++) 
					  {
					   printf("fv[%i]=%g\n",i,fv[i]);
					   printf("fb[%i]=%g\n",i,fb[i]);
					  }


                                         // ------------------------------------------------------------------------------	 		
                                         // end Todor Method Area


}


int main(void)
{
//point(double Qx, double Qy, double Qz, double vx, double vy, double vz, double bx, double by, double bz)
point(0.0, 0.0, 0.65, -1.0, 1.0, 0.0, 1.0, 1.0, 0.0);


return 0;
}
