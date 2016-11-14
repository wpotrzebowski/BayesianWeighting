#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>

struct mc_params { gsl_vector *saxs_ens;
                gsl_vector *saxs_exp;
                gsl_vector *err_saxs;
                double saxs_scale;
                gsl_vector *w_pre;
		        int k;
		        int N;};


double Energy(void *x, size_t k)
{
    //Data imports
    gsl_vector *saxs_ens = (gsl_vector *) (x->saxsEnsPtr);
    gsl_vector *saxs_exp = (gsl_vector *) (x->saxsExpPtr);
    gsl_vector *err_saxs = (gsl_vector *) (x->saxsErrPtr);
    gsl_vector *w_pre = (gsl_matrix *) (x->wPrePtr);
    double saxs_scale = x->saxsScale;

	double fit_prior = 1.0, fit_saxs = 1.0;

	for( int i = 0; i< N; i++) { fit_saxs *=
	exp( -(pow( saxs_scale*gsl_vector_get(saxs_ens,i) - gsl_vector_get(saxs_exp,i),2)/
	pow(gsl_vector_get(err_saxs,i),2))); }

	//for( int i = 0; i < k; i++) { fit_prior *= gsl_vector_get(w_pre,i); }
    //fit_prior *=  gsl_sf_lngamma(k/2)/k*gsl_sf_lngamma(0.5);
	return fit_saxs;
}

void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
  printf ("error  = % .6f = %.2g sigma\n", result - exact,
          fabs (result - exact) / error);
}

void mc_integrate(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs,
                double saxs_scale, gsl_vector *w_pre,
		        int k, int N) {
  double res, err;
  double *xl, *xu;

  xl = (double * ) malloc( k * sizeof( double ));
  xu = (double * ) malloc( k * sizeof( double ));

  for(int i = 0; i < k; i++) {
    xl[i] = 0.0;
    xu[i] = 1.0;
  }
  const gsl_rng_type *T;
  gsl_rng *r;

  struct mc_params params = {saxs_ens, saxs_exp, err_saxs, w_pre, saxs_scale, N};

  gsl_monte_function F;
  F.f = &Energy;
  F.k = k;
  F.params = &params;
  size_t calls = 500000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_monte_plain_state *s = gsl_monte_plain_alloc (k);
  gsl_monte_plain_integrate (&G, xl, xu, 3, calls, r, s,
                               &res, &err);
  gsl_monte_plain_free (s);

  display_results ("plain", res, err);

  /*
  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G, xl, xu, 3, calls, r, s,
                               &res, &err);
    gsl_monte_miser_free (s);

    display_results ("miser", res, err);
  }

  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
                               &res, &err);
    display_results ("vegas warm-up", res, err);

    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                                   &res, &err);
        printf ("result = % .6f sigma = % .6f "
                "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (s);
  }
  */
  gsl_rng_free (r);

}