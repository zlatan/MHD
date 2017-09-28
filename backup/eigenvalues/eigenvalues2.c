#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


int main (void)
{
  double w=-2.0/3;
  double Q;
  for(Q=0.0;Q<1.17;Q+=0.01)
  {
  double data[] = { 0.0, 0.0, -Q, 0.0,
                    1.0, 0.0, 0.0, -Q,
                     Q , 0.0, 0.0, 2*w,
		    0.0,  Q ,-2*w-1,0 };

  gsl_matrix_view m = gsl_matrix_view_array (data, 4, 4);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (4);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (4, 4);

  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (4);
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  
  {
    int i;
    double tmp[7];
    double vx[7],vy[7];
    for (i = 0; i < 4; i++)
      { 
        gsl_complex eval_i = gsl_vector_complex_get (eval, i);
	gsl_complex labdasq = gsl_complex_mul(eval_i,eval_i);
	tmp[i] = GSL_REAL(labdasq);
        gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
	
        vx[i]=GSL_REAL(gsl_vector_complex_get(&evec_i.vector, 0));
        vy[i]=GSL_REAL(gsl_vector_complex_get(&evec_i.vector, 2));

        //printf ("value -> %g\n", eval_i*eval_i);

        //printf ("eigenvector = \n");
        //gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
      }
    tmp[7]=tmp[0];
    for (i = 0; i < 4; i++)
      { 
       if ( tmp[i] > tmp[7] )
	{	
	tmp[7] = tmp[i];
	vx[7]  = vx[i];
        vy[7]  = vy[i];
        }
      }
     printf ("%g %g %g\n",Q,vx[7],vy[7]);
  }

  gsl_vector_free (eval);
  gsl_matrix_free (evec);
  }
  return 0;
}
