#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_statistics.h>
#include <time.h>
#include <omp.h>
#include "kde.h"

void progress_bar(float progress);

double ientropy(const gsl_vector *w, int k);

double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k);

void find_square_root(gsl_vector *w_ens, gsl_vector *w_ens1, double ct, double ct_prim, int k);

void ilr(int k, gsl_vector *in, gsl_vector *out);

void norm(gsl_vector *v);

void ilrInv( int k, gsl_vector *jerk, gsl_matrix *U, gsl_vector *in, gsl_vector *out);

double Force(gsl_vector *h_ens, gsl_vector *h_pre, int k);

double Energy(gsl_vector *h_ens, gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs,
	gsl_vector *cs_ens, gsl_vector *cs_exp, gsl_vector *cs_err_exp, gsl_vector *cs_err_pre,
	gsl_vector *h_pre, double saxs_scale,
	double f, int k, int N, int n, double T);

void Update(gsl_vector *h_ens,
	gsl_vector *w_ens,
	gsl_vector *saxs_ens, gsl_matrix *saxs_pre,
	gsl_vector *cs_ens, gsl_matrix *cs_pre,
	gsl_matrix *U, gsl_vector *jerk, int k);

void RandomStepH(gsl_vector *h_ens_current, gsl_vector *h_ens_trial,
	gsl_vector *w_ens_trial,
	gsl_vector *saxs_ens_trial, gsl_matrix *saxs_pre,
	gsl_vector *cs_ens_trial, gsl_matrix *cs_pre,
	gsl_matrix *U, gsl_vector *jerk, gsl_rng *r, double size, int k);

double SaxsScaleMean(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N);

double SaxsScaleStandardDeviation(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N, double T);
