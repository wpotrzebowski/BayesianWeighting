#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_poly.h>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include "gsl_siman.h"

typedef struct {
        size_t size;
        double *alphas;
        void *saxsExpPtr;
        void *saxsErrPtr;
        void *saxsEnsPtr;
        void *saxsPrePtr;
        void *saxsMixPtr;
        void *saxsScale;
        size_t numberProcs;
	void *OligomericSpecies;
        void *Concentration;
        double MonomerMass;
        size_t OligomerOrder;
	size_t numberOfCurves;
        } block;

block * block_alloc(size_t n);

void block_free(block * t);

void block_copy(void *inp, void *outp);

void * block_copy_construct(void *xp);

void block_destroy(void *xp);

double ientropy(const gsl_vector *w, int k);
double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k);
void polySolver (int order, double cmass_ratio, gsl_vector *oligomeric_states, gsl_matrix *kconsts,  double roots[]);
void find_poly_root(gsl_vector *w_ens, gsl_vector *w_ens_prim, double ct, double ct_prim,
        double monomerMass, int k, int order, gsl_vector *oligomeric_species );
double SaxsScaleMean(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N);
double SaxsScaleStandardDeviation(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N, double T);
double L_function(void *xp);
double L_distance(void *xp, void *yp);
void L_print (void *xp);
void L_take_step(const gsl_rng * r, void *xp, double step_size);
void run_vbw(const int &again, const int &k, const std::string &mdfile,
        const int &N, const std::string &presaxsfile,
        const int &Ncurves, const std::string &curvesfile,
        const std::string &filelist, const double &monomerMass,
        const std::string &outfile, const int &nprocs, const double &w_cut);

