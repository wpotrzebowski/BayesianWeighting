#include "VBW_sc.hh"


double ModelEvidenceEnergy(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs,
                double saxs_scale, int N)
{
	double fit_prior = 1.0, fit_saxs = 1.0;

	for( int i = 0; i< N; i++) { fit_saxs *=
	exp( -(pow( saxs_scale*gsl_vector_get(saxs_ens,i) - gsl_vector_get(saxs_exp,i),2)/
	pow(gsl_vector_get(err_saxs,i),2))); }

	return fit_saxs;
}


double mc_integrate(gsl_matrix *saxs_pre, gsl_vector *saxs_exp,
                    gsl_vector *err_saxs, int k, int N) {

  double energy_final;
  double *alphas;
  double *samples;
  double *alpha_ens;
  size_t Ntrials = 10000;

  alphas = (double * ) malloc( k * sizeof( double ));
  samples = (double * ) malloc( k * sizeof( double ));
  alpha_ens = (double * ) malloc( k * sizeof( double ));
  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_vector *weights = gsl_vector_alloc(k);
  gsl_vector *saxs_ens = gsl_vector_alloc(N);

  for (int i = 0; i<k; i++) alphas[i] = 0.5;

  double energy_trial=0.0;
  double saxs_scale = 0.0;
  double alpha_zero;
  for (int i=0; i<Ntrials; i++) {
    gsl_ran_dirichlet(r, k, alphas, samples);
    alpha_zero = 0.0;
    for (int j = 0; j < k; j++) {
        alpha_zero += samples[j];
	}

	for( int j = 0; j< N; j++) {
	    alpha_ens[j] = 0.0;
	    for (int l = 0; l < k; l++) {
		    alpha_ens[j]+=gsl_matrix_get(saxs_pre,j,l)*samples[l];
	    }
	    gsl_vector_set(weights, j, alpha_ens[j]/alpha_zero);
    }
    saxs_scale = SaxsScaleMean(weights,saxs_exp,err_saxs,N);
    gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, weights, 0.0, saxs_ens);
    energy_trial+=ModelEvidenceEnergy(saxs_ens,saxs_exp,err_saxs,saxs_scale,N);
  }
  energy_final/=Ntrials;
  gsl_rng_free (r);
  return energy_final;

}