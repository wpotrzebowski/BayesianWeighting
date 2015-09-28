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
#include <gsl/gsl_siman.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <sys/time.h>
#include <omp.h>

typedef struct {
        size_t size;
        double *alphas;
        void *saxsExpPtr;
        void *saxsErrPtr;
        void *saxsEnsPtr;
        void *saxsPrePtr;
        void *saxsMixPtr;
        double saxsScale;
        int numberProcs;
        } block;

block * block_alloc(size_t n);

void block_free(block * t);

void block_copy(void *inp, void *outp);

void * block_copy_construct(void *xp);

void block_destroy(void *xp);

double ientropy(const gsl_vector *w, int k);
double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k);
void find_square_root(gsl_vector *w_ens, gsl_vector *w_ens1, double ct, double ct_prim, int k);
double SaxsScaleMean(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N);
double SaxsScaleStandardDeviation(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N, double T);
double L_function(void *xp);
double L_distance(void *xp, void *yp);
void L_print (void *xp);
void L_take_step(const gsl_rng * r, void *xp, double step_size);
