#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

int
print_state (size_t iter, gsl_multiroot_fsolver * s, int n)
{
  
  printf ("iter = %3u",iter);
  for (int j=0; j<n; j++) {
	printf ("% .3f",gsl_vector_get (s->x, j));
  }
  printf ("\n");

 /* printf ("iter = %3u x = % .3f % .3f  % .3f "
          "f(x) = % .3e % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2),
          gsl_vector_get (s->f, 0), 
          gsl_vector_get (s->f, 1),
	  gsl_vector_get (s->f, 2)
	);*/
}


struct rparams
  {
    gsl_vector *fd;
    double R;
    int nweights;
  };

int
rosenbrock_f (const gsl_vector * x, void *params, 
              gsl_vector * f)
{
  double R = ((struct rparams *) params)->R;
  int nweights = ((struct rparams *) params)->nweights;
  gsl_vector * fd = ((struct rparams *) params)->fd;
  
  double weight_sum(0.0);
  for (int j=0; j<nweights; j++) {
	weight_sum+=gsl_vector_get(x, j);
  }
  double fm_prim_2 = (1 - weight_sum)*(1 - weight_sum);

  for (int j=0; j<nweights; j++) {
	gsl_vector_set(f, j, R * gsl_vector_get(fd,j) * fm_prim_2 - gsl_vector_get(x,j));
  }
  /*const double fd0 = gsl_vector_get (fd, 0);
  const double fd1 = gsl_vector_get (fd, 1);
  const double fd2 = gsl_vector_get (fd, 2);

  const double fdprim0 = gsl_vector_get (x, 0);
  const double fdprim1 = gsl_vector_get (x, 1);
  const double fdprim2 = gsl_vector_get (x, 2);


  const double y0 = R * fd0 * (1 - fdprim0 - fdprim1 - fdprim2)*(1 - fdprim0 - fdprim1 - fdprim2) - fdprim0;
  const double y1 = R * fd1 * (1 - fdprim0 - fdprim1 - fdprim2)*(1 - fdprim0 - fdprim1 - fdprim2) - fdprim1;
  const double y2 = R * fd2 * (1 - fdprim0 - fdprim1 - fdprim2)*(1 - fdprim0 - fdprim1 - fdprim2) - fdprim2;

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1); 
  gsl_vector_set (f, 2, y2);*/

  return GSL_SUCCESS;
}

int
main (void)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  double concentration1 = 2.3;
  double concentration2 = 4.6; 
  double fm = 0.4;
  double R = concentration2/(concentration1*fm*fm);
  int status;
  size_t i, iter = 0;

  const size_t n = 500;
  gsl_vector *fd = gsl_vector_alloc (n);
  
  for (int j=0; j<n; j++) {
	double step = j;
  	gsl_vector_set (fd, j, (j+1)/500.0);
  }
  struct rparams p = { fd, R, n };
  gsl_multiroot_function f = {&rosenbrock_f, n, &p};

  gsl_vector *x = gsl_vector_alloc (n);
  for (int j=0; j<n; j++) { 
  	gsl_vector_set (x, j, 0.001);
  }
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);

  print_state (iter, s, n);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      print_state (iter, s, n);

      if (status)   /* check if solver is stuck */
        break;

      status = 
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 10000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return 0;
}
