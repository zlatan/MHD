#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

int main (void)
{
  double Qa=0.6;
  double nx=0.2;
  double ny=0.4;
  double nz=0.3;
  double w=-3.0/2;
  double data[] = { 0.0, 0.0, 0.0, -Qa, 0.0, 0.0,
                    1.0, 0.0, 0.0, 0.0, -Qa, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, -Qa, 
                    Qa , 0.0, 0.0, 2*ny*nx*(w+1)        , -2*nx*nx*w+2*w , 0,
		    0.0, Qa , 0.0, 2*ny*ny*(w+1)-(2*w+1), -2*nx*ny*w     , 0,
		    0.0, 0.0, Qa , 2*ny*nz*(w+1)        , -2*nx*nz*w	 , 0 };

  gsl_matrix_view m = gsl_matrix_view_array (data, 6, 6);
  gsl_vector *eval = gsl_vector_alloc (6);
  gsl_matrix *evec = gsl_matrix_alloc (6, 6);
  gsl_eigen_symmv_workspace * ws = gsl_eigen_symmv_alloc (6);
  gsl_eigen_symmv (&m.matrix, eval, evec, ws);
  gsl_eigen_symmv_free (ws);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  
  {
    int i;

    for (i = 0; i < 6; i++)
      {
        double eval_i = gsl_vector_get (eval, i);
        gsl_vector_view evec_i = gsl_matrix_column (evec, i);
        printf ("eigenvalue = %g\n", eval_i);
        printf ("eigenvector = \n");
        gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
      }
  }

  gsl_vector_free (eval);
  gsl_matrix_free (evec);

  return 0;
}
