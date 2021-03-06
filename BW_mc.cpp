#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_statistics.h>
#include <omp.h>
#include "kde.h"

using namespace std;

// Constant defintions
const double pi = M_PI;

double ientropy(const gsl_vector *w, int k);
double ientropy(const gsl_vector *w, int k) {
	double ie = 0.0;
	for (int i=0; i<k; i++) 
		ie -= gsl_vector_get(w,i)*log2(gsl_vector_get(w,i));
	return ie;
}

double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k);
double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k) {

	double jsd;
	gsl_vector *w_c = gsl_vector_alloc(k);
	gsl_vector_memcpy(w_c, w_b);
	gsl_vector_add(w_c,w_a);
	gsl_vector_scale(w_c,0.5);
	jsd = ientropy(w_c,k) - 0.5*ientropy(w_a,k) - 0.5*ientropy(w_b,k);
	return jsd;
}


void find_square_root(gsl_vector *w_ens, gsl_vector *w_ens1, double ct, double ct_prim, int k);
void find_square_root(gsl_vector *w_ens, gsl_vector *w_ens1, double ct, double ct_prim, int k)
{
	double ksum;
	double w_prim_sum;
	double cn_inv, wm2_inv, cm_prim, cn_prim_inv, w_prim, kd;
	gsl_vector *kconst = gsl_vector_alloc(k-1);

	double wm = gsl_vector_get(w_ens,k-1);
	if (wm< 0.0001) wm = 0.0001;

	cn_inv = (2-wm)/ct;
	wm2_inv = 1/(wm*wm);
        	
	ksum = 0;
	for(int i = 0; i < (k-1); i++) {
		kd = gsl_vector_get(w_ens,i)*cn_inv*wm2_inv;
    		gsl_vector_set(kconst,i,kd);
    		ksum +=kd;
	}

	cm_prim = (sqrt(1+8*ksum*ct_prim)-1)/(4*ksum);
	cn_prim_inv = 2/(ct_prim+cm_prim);
	w_prim_sum = 0;
	for(int i = 0; i < (k-1); i++) {
		w_prim = gsl_vector_get(kconst,i)*cm_prim*cm_prim*cn_prim_inv;	
		gsl_vector_set(w_ens1,i,w_prim);
		w_prim_sum +=w_prim;
        }
	gsl_vector_set(w_ens1,k-1,cm_prim*cn_prim_inv);
	gsl_vector_free(kconst);
}

// Function defintions
//Equation 4
//Setting h vector for isomorphic transformation from Simplex to Euclidian space
void ilr(int k, gsl_vector *in, gsl_vector *out);
void ilr(int k, gsl_vector *in, gsl_vector *out)
{
	double j = 0.0, temp = 0.0;
	for(int i = 0; i < (k-1); i++) 
	{
		j = i +1.0; temp = 0.0;
		for(int l = 0; l < j; l++) { temp += log(gsl_vector_get(in,l)); } 
		gsl_vector_set(out, i, (1.0/sqrt( j*(j+1) ))*(temp - j*log(gsl_vector_get(in,j))) );		
	}
}

void norm(gsl_vector *v);
void norm(gsl_vector *v) { gsl_vector_scale(v,1.0/gsl_blas_dasum(v)); } //1 over absolute sum

void ilrInv( int k, gsl_vector *jerk, gsl_matrix *U, gsl_vector *in, gsl_vector *out); 
void ilrInv( int k, gsl_vector *jerk, gsl_matrix *U, gsl_vector *in, gsl_vector *out)
{
	gsl_blas_dgemv(CblasNoTrans,1.0,U,in,0.0,jerk);
	for (int i = 0; i < k; i++)  { gsl_vector_set(out,i,exp(gsl_vector_get(jerk,i))); }
	norm(out);
}

//Equation 5. Calculates prior based on h vectors difference
double Force(gsl_vector *h_ens, gsl_vector *h_pre, int k); 
double Force(gsl_vector *h_ens, gsl_vector *h_pre, int k)
{
	double fit_prior = 0.0;
	for( int i = 0; i < k-1; i++) { fit_prior += pow( gsl_vector_get(h_ens,i) - gsl_vector_get(h_pre,i), 2); }
	return 0.5 * fit_prior;
}

double Energy(gsl_vector *h_ens, gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, 
		gsl_vector *saxs_ens1, gsl_vector *saxs_exp1, gsl_vector *err_saxs1, 
		gsl_vector *saxs_ens2, gsl_vector *saxs_exp2, gsl_vector *err_saxs2, 
		gsl_vector *h_pre, double saxs_scale, double saxs_scale1, double saxs_scale2,
		double f, int k, int N, int N1, int N2, double T);

double Energy(gsl_vector *h_ens, gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, 
		gsl_vector *saxs_ens1, gsl_vector *saxs_exp1, gsl_vector *err_saxs1, 
		gsl_vector *saxs_ens2, gsl_vector *saxs_exp2, gsl_vector *err_saxs2, 
		gsl_vector *h_pre, double saxs_scale, double saxs_scale1, double saxs_scale2,
		double f, int k, int N, int N1, int N2,  double T)
{
	double fit_prior = 0.0, fit_saxs = 0.0, fit_saxs1 = 0.0, fit_saxs2 = 0.0;
	for( int i = 0; i< N; i++) { fit_saxs += (pow( saxs_scale*gsl_vector_get(saxs_ens,i) - gsl_vector_get(saxs_exp,i), 2) / gsl_vector_get(err_saxs,i) ); }
	for( int i = 0; i< N1; i++) { fit_saxs1 += (pow( saxs_scale1*gsl_vector_get(saxs_ens1,i) - gsl_vector_get(saxs_exp1,i), 2) / gsl_vector_get(err_saxs1,i) ); }
	for( int i = 0; i< N2; i++) { fit_saxs2 += (pow( saxs_scale2*gsl_vector_get(saxs_ens2,i) - gsl_vector_get(saxs_exp2,i), 2) / gsl_vector_get(err_saxs2,i) ); }
	for( int i = 0; i < k-1; i++) { fit_prior += pow( gsl_vector_get(h_ens,i) - gsl_vector_get(h_pre,i), 2) * f; }
	return 0.5*(fit_saxs + fit_saxs1 + fit_saxs2 + fit_prior)/T;
	//return 0.5*( fit_prior)/T;
	//1/T comes from the multiple replica exchnages
}

void Update(gsl_vector *h_ens, 
	gsl_vector *w_ens,  
	gsl_vector *w_ens1,
	gsl_vector *w_ens2, 
	gsl_vector *saxs_ens, gsl_matrix *saxs_pre, 
	gsl_vector *saxs_ens1, gsl_matrix *saxs_pre1, 
	gsl_vector *saxs_ens2, gsl_matrix *saxs_pre2, 
	gsl_matrix *U, gsl_vector *jerk, int k);
void Update(gsl_vector *h_ens, 
	gsl_vector *w_ens, 
	gsl_vector *w_ens1, 
	gsl_vector *w_ens2, 
	gsl_vector *saxs_ens, gsl_matrix *saxs_pre, 
	gsl_vector *saxs_ens1, gsl_matrix *saxs_pre1, 
	gsl_vector *saxs_ens2, gsl_matrix *saxs_pre2, 
	gsl_matrix *U, gsl_vector *jerk, int k)
{
	double cm, cm_prim;
	ilrInv(k,jerk,U,h_ens,w_ens);
	//First concentartion transition
	cm = 4.6;
        cm_prim = 2.3;
	find_square_root(w_ens,w_ens1,cm,cm_prim,k);
	//KAUpdate(w_ens,w_ens1,w_pre,w_pre1,cm,cm_prim,k);
        cm = 4.6;
        cm_prim = 1.15;
	find_square_root(w_ens,w_ens2,cm,cm_prim,k);	
	//KAUpdate(w_ens,w_ens2,w_pre,w_pre2,cm,cm_prim,k);

	//These functions compute the matrix-vector product and sum // y = 1.0*A*x
	gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_ens, 0.0, saxs_ens);
	gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre1, w_ens1, 0.0, saxs_ens1);
	gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre2, w_ens2, 0.0, saxs_ens2);
}

//This function most likely can be run independently of weights
void RandomStepH(gsl_vector *h_ens_current, gsl_vector *h_ens_trial,
	gsl_vector *w_ens_trial,  
	gsl_vector *w_ens_trial1, 
	gsl_vector *w_ens_trial2, 
	gsl_vector *saxs_ens_trial, gsl_matrix *saxs_pre, 
	gsl_vector *saxs_ens_trial1, gsl_matrix *saxs_pre1, 
	gsl_vector *saxs_ens_trial2, gsl_matrix *saxs_pre2, 
	gsl_matrix *U, gsl_vector *jerk, gsl_rng *r, double size, int k);
void RandomStepH(gsl_vector *h_ens_current, gsl_vector *h_ens_trial, 
	gsl_vector *w_ens_trial, 
	gsl_vector *w_ens_trial1,
	gsl_vector *w_ens_trial2,
	gsl_vector *saxs_ens_trial, gsl_matrix *saxs_pre, 
	gsl_vector *saxs_ens_trial1, gsl_matrix *saxs_pre1, 
	gsl_vector *saxs_ens_trial2, gsl_matrix *saxs_pre2, 
	gsl_matrix *U, gsl_vector *jerk, gsl_rng *r, double size, int k)
{
	for(int i = 0; i < k-1; i++) { 
		gsl_vector_set(h_ens_trial, i, gsl_vector_get(h_ens_current,i) + gsl_ran_gaussian(r,size) ); 
	}
	Update(h_ens_trial, 
		w_ens_trial, 
		w_ens_trial1,  
		w_ens_trial2, 
		saxs_ens_trial,saxs_pre,
		saxs_ens_trial1,saxs_pre1,
		saxs_ens_trial2,saxs_pre2,
		U,jerk,k);
}

double SaxsScaleMean(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N); 
double SaxsScaleMean(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N)
{
	double tempa = 0.0, tempb = 0.0;
	for( int i = 0; i< N; i++) {
		tempa += gsl_vector_get(saxs_ens,i)*gsl_vector_get(saxs_exp,i)/gsl_vector_get(err_saxs,i);
		tempb += pow(gsl_vector_get(saxs_ens,i),2.0)/gsl_vector_get(err_saxs,i);
	}
	return tempa/tempb;
}

double SaxsScaleStandardDeviation(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N, double T); 
double SaxsScaleStandardDeviation(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N, double T)
{
	double temp = 0.0;
	for( int i = 0; i< N; i++) { 
		temp += pow(gsl_vector_get(saxs_ens,i),2.0)/gsl_vector_get(err_saxs,i); 
	}
	return sqrt(T/temp);
}

int main()
{
	int n,k,steps,equilibration,np,samples,skip,swap_frequency = 100, num_swaps, rep = 0, again = 0;
	char mdfile[80], outfile[80];
	//Number of measurements in first curve
	int N; 
	char presaxsfile[80], saxsfile[80], saxserrfile[80]; 
	//Number of measuremnets in second curve
	int N1;
        char presaxsfile1[80], saxsfile1[80], saxserrfile1[80];
	int N2;
        char presaxsfile2[80], saxsfile2[80], saxserrfile2[80];
	//Restart from already precaclculated vaules
	fscanf(stdin, "%d", &again); 
	//Number of processors/temperatures
	fscanf(stdin, "%d", &np); 
	//Number of strcutures in ensemble
	fscanf(stdin, "%d", &k); 
	//Prior weights
	fscanf(stdin, "%s", &mdfile);
	cout<<"Reading prior values"<<std::endl; 
	////////////////////////////
	cout<<"Reading 1st scatteirng curve"<<std::endl;
	//Number of SAXS measurements in curve 1
	fscanf(stdin, "%d", &N); 
	fscanf(stdin, "%s", &presaxsfile); 
	fscanf(stdin, "%s", &saxsfile); 
	fscanf(stdin, "%s", &saxserrfile); 

	cout<<"Reading 2nd scatteirng curve"<<std::endl;	
	//Number of SAXS measurements in curve 2
        fscanf(stdin, "%d", &N1);
        fscanf(stdin, "%s", &presaxsfile1);
        fscanf(stdin, "%s", &saxsfile1);
        fscanf(stdin, "%s", &saxserrfile1);

	cout<<"Reading 3rd scatteirng curve"<<std::endl;
	//Number of SAXS measurements in curve 3
        fscanf(stdin, "%d", &N2);
        fscanf(stdin, "%s", &presaxsfile2);
        fscanf(stdin, "%s", &saxsfile2);
        fscanf(stdin, "%s", &saxserrfile2);

	//Running params
	fscanf(stdin, "%s", &outfile); 
	fscanf(stdin, "%d", &equilibration); 
	fscanf(stdin, "%d", &steps); 
	fscanf(stdin, "%d", &samples); 
		
	double saxs_scale_current[np], saxs_scale_current1[np], saxs_scale_current2[np];
	double h[k], f[np],  accepted[np], step_size[np], temperature[np], swaps_accepted[np];
	double energy_current[np], energy_trial[np]; 	
	double fL = 0.001, dh = 0.001, temp = 0.0, j = 0.0;

	gsl_matrix *tostart = gsl_matrix_alloc(np, k+2), 
		*U = gsl_matrix_alloc(k,k-1), 
		*saxs_pre = gsl_matrix_alloc(N,k), 
		*saxs_pre1 = gsl_matrix_alloc(N1,k),
		*saxs_pre2 = gsl_matrix_alloc(N2,k),
		//Two vectors of weights plus some sampling info
		*memory = gsl_matrix_alloc(samples,3*k+5),
		*basis = gsl_matrix_alloc(k-1,k),
		*weight_samples = gsl_matrix_alloc(steps,k),
		*weight_samples1 = gsl_matrix_alloc(steps,k),
		*weight_samples2 = gsl_matrix_alloc(steps,k);

	gsl_vector *saxs_exp = gsl_vector_alloc(N),
		*err_saxs = gsl_vector_alloc(N),
		*saxs_exp1 = gsl_vector_alloc(N1),
                *err_saxs1 = gsl_vector_alloc(N1),
		*saxs_exp2 = gsl_vector_alloc(N2),
                *err_saxs2 = gsl_vector_alloc(N2),
		*jerk[np],
		*w_pre = gsl_vector_alloc(k),
		*h_pre = gsl_vector_alloc(k-1),
		*h_mid = gsl_vector_alloc(k-1),
		*w_ens_current[np],
		*w_ens_current1[np],
		*w_ens_current2[np],
		*h_ens_current[np],
		*saxs_ens_current[np],
		*saxs_ens_current1[np],
		*saxs_ens_current2[np],
		*pdm_vec[np],
		*pdm_vec1[np],
		*pdm_vec2[np],
		*w_ens_trial[np],
		*w_ens_trial1[np],
		*w_ens_trial2[np],
		*h_ens_trial[np],
		*saxs_ens_trial[np],
		*saxs_ens_trial1[np],
		*saxs_ens_trial2[np],
		*bayesian_weight1 = gsl_vector_alloc(k),
		*bayesian_weight2 = gsl_vector_alloc(k),
		*bayesian_weight3 = gsl_vector_alloc(k),
		*bayesian_weight1_current = gsl_vector_alloc(k),
                *bayesian_weight2_current = gsl_vector_alloc(k),
                *bayesian_weight3_current = gsl_vector_alloc(k),
		*w_ens_last_accepted = gsl_vector_alloc(k),
		*w_ens_last_accepted1 = gsl_vector_alloc(k),
		*w_ens_last_accepted2 = gsl_vector_alloc(k);


	gsl_vector_set_zero(bayesian_weight1);
	gsl_vector_set_zero(bayesian_weight2);
	gsl_vector_set_zero(bayesian_weight3);		 
	// intialize variables //
	skip = steps / samples;
	num_swaps = steps / swap_frequency;
	omp_set_num_threads(np); 

	for(int i = 0; i < np; i++) 
	{ 
		if (np == 1) temperature[i] = 1;
		else temperature[i] = pow( 1.5,float(i)/(float(np) - 1.0) );
		cout<<"Temperature: "<<temperature[i]<<std::endl; 
		accepted[i] = 0.0; 
		step_size[i] = 0.1; 
		swaps_accepted[i] = 0.0;
		jerk[i] = gsl_vector_alloc(k);
		w_ens_current[i] = gsl_vector_alloc(k); 
		w_ens_current1[i] = gsl_vector_alloc(k);
		w_ens_current2[i] = gsl_vector_alloc(k);
		h_ens_current[i] = gsl_vector_alloc(k-1); 
		saxs_ens_current[i] = gsl_vector_alloc(N);
		saxs_ens_current1[i] = gsl_vector_alloc(N1);
		saxs_ens_current2[i] = gsl_vector_alloc(N2);
		pdm_vec[i] = gsl_vector_alloc(N);
		pdm_vec1[i] = gsl_vector_alloc(N1);
		pdm_vec2[i] = gsl_vector_alloc(N2);
		w_ens_trial[i] = gsl_vector_alloc(k); 
		w_ens_trial1[i] = gsl_vector_alloc(k);
		w_ens_trial2[i] = gsl_vector_alloc(k);
		h_ens_trial[i] = gsl_vector_alloc(k-1); 
		saxs_ens_trial[i] = gsl_vector_alloc(N);
		saxs_ens_trial1[i] = gsl_vector_alloc(N1);
		saxs_ens_trial2[i] = gsl_vector_alloc(N2);
	}
	// Read in data from files //
	FILE * inFile = fopen(presaxsfile,"r"); gsl_matrix_fscanf(inFile,saxs_pre);fclose(inFile);
	inFile = fopen(saxsfile,"r"); gsl_vector_fscanf(inFile,saxs_exp); fclose(inFile);
	inFile = fopen(saxserrfile,"r"); gsl_vector_fscanf(inFile,err_saxs); fclose(inFile);
	inFile = fopen(presaxsfile1,"r"); gsl_matrix_fscanf(inFile,saxs_pre1); fclose(inFile);
        inFile = fopen(saxsfile1,"r"); gsl_vector_fscanf(inFile,saxs_exp1); fclose(inFile);
        inFile = fopen(saxserrfile1,"r"); gsl_vector_fscanf(inFile,err_saxs1); fclose(inFile);
	inFile = fopen(presaxsfile2,"r"); gsl_matrix_fscanf(inFile,saxs_pre2); fclose(inFile);
        inFile = fopen(saxsfile2,"r"); gsl_vector_fscanf(inFile,saxs_exp2); fclose(inFile);
        inFile = fopen(saxserrfile2,"r"); gsl_vector_fscanf(inFile,err_saxs2); fclose(inFile);
	inFile = fopen(mdfile,"r"); gsl_vector_fscanf(inFile,w_pre); fclose(inFile);
	if(again == 1){ inFile = fopen("restart.dat","r"); gsl_matrix_fscanf(inFile,tostart); fclose(inFile); }

	//Create matrix for basis transformation//
	//First part of equation (4)//
	temp = 0.0, j = 0.0;
	for(int q = 0; q < (k-1); q++) {
		for(int i = 0; i < k; i++) {
			j = q + 1.0;
			if(i<j) { temp = 1/j; }
			else if(i == j) { temp = -1.0; }
			else { temp = 0.0; }
			gsl_matrix_set(basis,q,i,sqrt( j/(j+1)) * temp);
		}
	}

	for(int i = 0; i < k-1; i++) { 
		for(int j = 0 ; j < k; j++) { 
			gsl_matrix_set(U,j,i,gsl_matrix_get(basis,i,j)); 
		} 
	}	
	// initialize random number generators //
	const gsl_rng_type *K; 
	gsl_rng *r[np]; 
	gsl_rng_env_setup(); 
	K = gsl_rng_default;
	
	for(int i = 0; i < np; i++) { 
		r[i] = gsl_rng_alloc(K); 
		gsl_rng_set(r[i],time(NULL)+i); 
	}
	
	// initialize variables //
	//Setting h-vector//
	ilr(k,w_pre,h_pre);
	gsl_vector_set_zero(h_mid);
	for(rep = 0; rep < np; rep++)
	{
		gsl_vector_memcpy(h_ens_current[rep],h_mid);
		//This function returns a random variate from the exponential distribution with mean mu
		f[rep] = gsl_ran_exponential(r[rep], 1.0/Force(h_ens_current[rep],h_pre,k)) + fL;
		saxs_scale_current[rep] = 1.0;
		saxs_scale_current1[rep] = 1.0;
		saxs_scale_current2[rep] = 1.0;
		RandomStepH(h_ens_current[rep],h_ens_current[rep],
			w_ens_current[rep], 
			w_ens_current1[rep],
			w_ens_current2[rep],
			saxs_ens_current[rep],saxs_pre,
			saxs_ens_current1[rep],saxs_pre1,
			saxs_ens_current2[rep],saxs_pre2,
			U,jerk[rep],r[rep],1.0,k);
	}

	cout<<"Variables have been initialized"<<std::endl;	
	if(again==1)
	{
		for(int i = 0; i < np; i++)
		{
			for(int j = 0; j < k-1; j++) { 
				gsl_vector_set(h_ens_current[i],j,gsl_matrix_get(tostart,i,j)); 
			}
			f[i] = gsl_matrix_get(tostart,i,k-1); 
			saxs_scale_current[i] = gsl_matrix_get(tostart,i,k); 
			saxs_scale_current1[i] = gsl_matrix_get(tostart,i,k);
			saxs_scale_current2[i] = gsl_matrix_get(tostart,i,k);
			step_size[i] = gsl_matrix_get(tostart,i,k+1);
			Update(h_ens_current[i],
				w_ens_current[i],
				w_ens_current1[i], 
				w_ens_current2[i], 
				saxs_ens_current[i],saxs_pre,
				saxs_ens_current1[i],saxs_pre1,
				saxs_ens_current2[i],saxs_pre2,U,jerk[i],k);
		}
	}
	
	//In general generates random variates with the given distribution and calculates energy based on these
	if(again != 1)
	{
		cout << "Equilibration" << endl;
		rep = 0;
		#pragma omp parallel for private(rep)
		for(rep = 0; rep < np; rep++)
		{
			for(int j = 0; j < equilibration; j++)
			{
				//Sampling over entire weights domain
				f[rep] = gsl_ran_exponential(r[rep], temperature[rep]/Force(h_ens_current[rep],h_pre,k)) + fL;

				//Sampling scaling factor for each scattering curve	
				saxs_scale_current[rep] = SaxsScaleMean(saxs_ens_current[rep],saxs_exp,err_saxs,N) + gsl_ran_gaussian(r[rep],SaxsScaleStandardDeviation(saxs_ens_current[rep],saxs_exp,err_saxs,N,temperature[rep]));
				saxs_scale_current1[rep] = SaxsScaleMean(saxs_ens_current1[rep],saxs_exp1,err_saxs1,N1) + gsl_ran_gaussian(r[rep],SaxsScaleStandardDeviation(saxs_ens_current1[rep],saxs_exp1,err_saxs1,N1,temperature[rep]));
				saxs_scale_current2[rep] = SaxsScaleMean(saxs_ens_current2[rep],saxs_exp2,err_saxs2,N2) + gsl_ran_gaussian(r[rep],SaxsScaleStandardDeviation(saxs_ens_current2[rep],saxs_exp2,err_saxs2,N2,temperature[rep]));

				RandomStepH(h_ens_current[rep],h_ens_trial[rep],
					w_ens_trial[rep],  
					w_ens_trial1[rep],
					w_ens_trial2[rep],
					saxs_ens_trial[rep],saxs_pre,
					saxs_ens_trial1[rep],saxs_pre1,
					saxs_ens_trial2[rep],saxs_pre2,
					U,jerk[rep],r[rep],step_size[rep],k);
				

				energy_current[rep] = Energy(h_ens_current[rep],
					saxs_ens_current[rep],saxs_exp,err_saxs,
					saxs_ens_current1[rep],saxs_exp1,err_saxs1,
					saxs_ens_current2[rep],saxs_exp2,err_saxs2,
					h_pre,saxs_scale_current[rep],saxs_scale_current1[rep],saxs_scale_current2[rep],
					f[rep],k,N,N1,N2,temperature[rep]);
		

				energy_trial[rep] = Energy(h_ens_trial[rep],
					saxs_ens_trial[rep],saxs_exp,err_saxs,
					saxs_ens_trial1[rep],saxs_exp1,err_saxs1,
					saxs_ens_trial2[rep],saxs_exp2,err_saxs2,
					h_pre,saxs_scale_current[rep],saxs_scale_current1[rep],saxs_scale_current2[rep],
					f[rep],k,N,N1,N2,temperature[rep]);
	
				//Monte Carlo accpeptance in terms of energies
				if(gsl_rng_uniform(r[rep]) <= exp(-energy_trial[rep] + energy_current[rep]) )
				{
					gsl_vector_memcpy(h_ens_current[rep],h_ens_trial[rep]);
					Update(h_ens_current[rep],
						w_ens_current[rep],
						w_ens_current1[rep],
						w_ens_current2[rep],
						saxs_ens_current[rep],saxs_pre,
						saxs_ens_current1[rep],saxs_pre1,
						saxs_ens_current2[rep],saxs_pre2,
						U,jerk[rep],k);
	
					accepted[rep] += 1.0;
					energy_current[rep] = energy_trial[rep];
				}	
				//Proposal distribution is tuned so that 25% of steps are accepted	
				if((j+1)%100 == 0) 
				{
					step_size[rep] = step_size[rep] * pow((1 - (0.24 - accepted[rep]/100.0)),2); 
					accepted[rep] = 0.0;
				}
			}
		}	
	}
	
	for( int i = 0; i < np; i++) { accepted[i] = 0.0; }
	cout << "Sampling" << endl;
	double jsd1_sum = 0.0;
	double jsd2_sum = 0.0;
	double jsd3_sum = 0.0;
	int sampling_step = 0;
	for(int z = 0; z < num_swaps; z++)
	{
		rep = 0;
		#pragma omp parallel for private(rep)
		for(rep = 0; rep < np; rep++)
		{
			for(int j = 0; j < swap_frequency; j++)
			{
				f[rep] = gsl_ran_exponential(r[rep], temperature[rep]/Force(h_ens_current[rep],h_pre,k)) + fL;
				
				saxs_scale_current[rep] = SaxsScaleMean(saxs_ens_current[rep],saxs_exp,err_saxs,N) + gsl_ran_gaussian(r[rep],SaxsScaleStandardDeviation(saxs_ens_current[rep],saxs_exp,err_saxs,N,temperature[rep]));
				saxs_scale_current1[rep] = SaxsScaleMean(saxs_ens_current1[rep],saxs_exp1,err_saxs1,N1) + gsl_ran_gaussian(r[rep],SaxsScaleStandardDeviation(saxs_ens_current1[rep],saxs_exp1,err_saxs1,N1,temperature[rep]));
				saxs_scale_current2[rep] = SaxsScaleMean(saxs_ens_current2[rep],saxs_exp2,err_saxs2,N2) + gsl_ran_gaussian(r[rep],SaxsScaleStandardDeviation(saxs_ens_current2[rep],saxs_exp2,err_saxs2,N2,temperature[rep]));

				RandomStepH(h_ens_current[rep],h_ens_trial[rep],
						w_ens_trial[rep], 
						w_ens_trial1[rep],
						w_ens_trial2[rep],
						saxs_ens_trial[rep],saxs_pre,
						saxs_ens_trial1[rep],saxs_pre1,
						saxs_ens_trial2[rep],saxs_pre2,
						U,jerk[rep],r[rep],step_size[rep],k);
				
				energy_current[rep] = Energy(h_ens_current[rep],saxs_ens_current[rep],saxs_exp,err_saxs,
						saxs_ens_current1[rep],saxs_exp1,err_saxs1,
						saxs_ens_current2[rep],saxs_exp2,err_saxs2,
						h_pre,saxs_scale_current[rep],saxs_scale_current1[rep],saxs_scale_current2[rep],
						f[rep],k,N,N1,N2,temperature[rep]);

				energy_trial[rep] = Energy(h_ens_trial[rep],saxs_ens_trial[rep],saxs_exp,err_saxs,
						saxs_ens_trial1[rep],saxs_exp1,err_saxs1,
						saxs_ens_trial2[rep],saxs_exp2,err_saxs2,
						h_pre,saxs_scale_current[rep],saxs_scale_current1[rep],saxs_scale_current2[rep],
						f[rep],k,N,N1,N2,temperature[rep]);

				if(gsl_rng_uniform(r[rep]) <= exp(-energy_trial[rep] + energy_current[rep]) )
				{	
					gsl_vector_memcpy(h_ens_current[rep],h_ens_trial[rep]);
					Update(h_ens_current[rep],
						w_ens_current[rep],
						w_ens_current1[rep],
						w_ens_current2[rep],
						saxs_ens_current[rep],saxs_pre,
						saxs_ens_current1[rep],saxs_pre1,
						saxs_ens_current2[rep],saxs_pre2,
						U,jerk[rep],k);
					accepted[rep] += 1.0;
					energy_current[rep] = energy_trial[rep];
					if(rep ==0) {
                                                for (int jind=0; jind<k; jind++) {
                                                        gsl_matrix_set(weight_samples,sampling_step,jind,gsl_vector_get(w_ens_current[0],jind));
							gsl_matrix_set(weight_samples1,sampling_step,jind,gsl_vector_get(w_ens_current1[0],jind));
							gsl_matrix_set(weight_samples2,sampling_step,jind,gsl_vector_get(w_ens_current2[0],jind));
                                                }
                                                gsl_vector_memcpy(w_ens_last_accepted,w_ens_current[0]);
						gsl_vector_memcpy(w_ens_last_accepted1,w_ens_current1[0]);
						gsl_vector_memcpy(w_ens_last_accepted2,w_ens_current2[0]);
                                                sampling_step++;
                                        }
				}	
			
				
				if(rep ==0)
				{
					int foo= z*swap_frequency + j + 1;
					double niter = 1.0/double(foo);
					double jsd1 = 0.0; 
					double jsd2 = 0.0;
					double jsd3 = 0.0;
					gsl_vector_add(bayesian_weight1,w_ens_current[0]);
					gsl_vector_add(bayesian_weight2,w_ens_current1[0]);
					gsl_vector_add(bayesian_weight3,w_ens_current2[0]);
					gsl_vector_memcpy(bayesian_weight1_current,bayesian_weight1);
					gsl_vector_memcpy(bayesian_weight2_current,bayesian_weight2);
					gsl_vector_memcpy(bayesian_weight3_current,bayesian_weight3);
					gsl_vector_scale(bayesian_weight1_current,niter);
					gsl_vector_scale(bayesian_weight2_current,niter);
					gsl_vector_scale(bayesian_weight3_current,niter);
					jsd1 = jensen_shannon_div(bayesian_weight1_current,w_ens_current[0],k);
					jsd2 = jensen_shannon_div(bayesian_weight2_current,w_ens_current1[0],k);
					jsd3 = jensen_shannon_div(bayesian_weight3_current,w_ens_current2[0],k);
					jsd1_sum += jsd1;
					jsd2_sum += jsd2;
					jsd3_sum += jsd3;
					if( foo % skip == 0)
					{
						cout<<"PED1 "<<sqrt(jsd1_sum*niter)<<std::endl;
						cout<<"PED2 "<<sqrt(jsd2_sum*niter)<<std::endl;
						cout<<"PED3 "<<sqrt(jsd3_sum*niter)<<std::endl;
						for( int l = 0; l < k; l++) { 
							gsl_matrix_set(memory,foo / skip -1, l , gsl_vector_get(bayesian_weight1_current,l));
							//gsl_matrix_set(memory,foo / skip -1, l , gsl_vector_get(w_ens_current[0],l)); 
						}

						for( int l = k; l < 2*k; l++) {
                                                        gsl_matrix_set(memory,foo / skip -1, l , gsl_vector_get(bayesian_weight2_current,l-k));
							//gsl_matrix_set(memory,foo / skip -1, l , gsl_vector_get(w_ens_current1[0],l-k));
                                                }
						for( int l = 2*k; l < 3*k; l++) {
                                                        gsl_matrix_set(memory,foo / skip -1, l , gsl_vector_get(bayesian_weight3_current,l-2*k));
							//gsl_matrix_set(memory,foo / skip -1, l , gsl_vector_get(w_ens_current2[0],l-2*k));
                                                }
						gsl_matrix_set(memory, foo/skip -1, 3*k, f[0]); 
						gsl_matrix_set(memory, foo/skip -1, 3*k+1, saxs_scale_current[0]);
						gsl_matrix_set(memory, foo/skip -1, 3*k+2, saxs_scale_current1[0]);
						gsl_matrix_set(memory, foo/skip -1, 3*k+3, saxs_scale_current2[0]); 
						gsl_matrix_set(memory, foo/skip -1, 3*k+4, energy_current[0]);
					}
				}
			}
		}
		// temperature swap //
		for( int i = 0; i < np; i++)
		{
			//Even-odd rule 
			if(z%2 == 0) { if(i%2 == 0 || i >= np-1) { continue; } }
			if(z%2 != 0) { if(i%2 != 0 || i >= np-1) { continue; } }
			int swapa = i, swapb = i+1;
			energy_current[swapa] = energy_current[swapa]*temperature[swapa]; energy_current[swapb] = energy_current[swapb]*temperature[swapb];
			if(gsl_rng_uniform(r[0]) < exp( ( (1.0/temperature[swapb]) - (1.0/temperature[swapa]) )*( energy_current[swapb] - energy_current[swapa])) )
			{
				gsl_vector_swap(h_ens_current[swapa],h_ens_current[swapb]); 

				gsl_vector_swap(w_ens_current[swapa],w_ens_current[swapb]); 
                                gsl_vector_swap(w_ens_current1[swapa],w_ens_current1[swapb]);
				gsl_vector_swap(w_ens_current2[swapa],w_ens_current2[swapb]);				

				gsl_vector_swap(saxs_ens_current[swapa],saxs_ens_current[swapb]);
				gsl_vector_swap(saxs_ens_current1[swapa],saxs_ens_current1[swapb]);
				gsl_vector_swap(saxs_ens_current2[swapa],saxs_ens_current2[swapb]);
				
				//TODO: Check if it doesn't make redundant thing
				temp = saxs_scale_current[swapa]; 
				saxs_scale_current[swapa] = saxs_scale_current[swapb]; 
				saxs_scale_current[swapb] = temp;
				
				temp = saxs_scale_current1[swapa];
                                saxs_scale_current1[swapa] = saxs_scale_current1[swapb];
                                saxs_scale_current1[swapb] = temp;

				temp = saxs_scale_current2[swapa];
                                saxs_scale_current2[swapa] = saxs_scale_current2[swapb];
                                saxs_scale_current2[swapb] = temp;

				temp = f[swapa]; 
				f[swapa] = f[swapb]; 
				f[swapb] = temp;
				swaps_accepted[i] += 1.0;
			}
		}
	}
	
	 //MDOEL SELECTION       
        cout<<"Starting Model Evidence with matrix of size: "<<sampling_step<<std::endl;
        long double logBMS;
	long double logKDE1;
	long double logKDE2;
	long double logKDE3;
        logKDE1 = log(kerneldensity(weight_samples, w_ens_last_accepted, sampling_step, k));
	logKDE2 = log(kerneldensity(weight_samples1, w_ens_last_accepted1, sampling_step, k));
	logKDE3 = log(kerneldensity(weight_samples2, w_ens_last_accepted2, sampling_step, k));
        //logKDE = log(kerneldensity(weight_samples,w_ens_last_accepted,sampling_step,k));
        logBMS = -energy_current[0];
        cout<<"Number of steps, PDM, KDE "<<logBMS<<" "<<logKDE1<<" :  "<<logKDE2<<" : "<<logKDE3<<std::endl;
        ofstream pdm("pdm.dat");
        pdm<<logBMS - logKDE1 - logKDE2 - logKDE3<<std::endl;
        pdm.close();

	// output //
	ofstream output(outfile);
	for (int i = 0; i < samples; i++) { for( int j = 0; j < 3*k+5; j++) { output << gsl_matrix_get(memory,i,j) << " "; if (j == 3*k+4) { output << endl; } } }
	output.close();

	ofstream restart("restart.dat");
	for(int i = 0; i < np; i++) 
	{ 
		for(int j = 0; j < k-1; j++) { restart << gsl_vector_get(h_ens_current[i],j) << " "; }
		restart << f[i] << " " << saxs_scale_current[i] << " " << saxs_scale_current1[i]<<" "<<saxs_scale_current2[i] << " " << step_size[i] << endl;
	}
	restart.close();
	
	for(int i = 0; i < np; i++)
	{
		cout << "chain: " << i << " temperature: " << temperature[i] << " percent steps accepted: " <<  accepted[i]/steps *100 << endl; 
		cout << step_size[i] << endl;
		cout << 2.0* swaps_accepted[i] / float(num_swaps) << endl << endl;
	}
	
	return 0;
}
