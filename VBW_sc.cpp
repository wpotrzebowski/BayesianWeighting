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
#include <gsl/gsl_siman.h>
#include <gsl/gsl_math.h>
#include <omp.h>

using namespace std;
const double pi = M_PI;

///////////////Structure which is needed for simulated annealing//////////// 
typedef struct {
        size_t size;
        double *alphas;
        void *saxsExpPtr;
	void *saxsErrPtr;
	void *saxsEnsPtr;
	void *saxsPrePtr;
	double saxsScale;
        } block;

block * block_alloc(size_t n) {
        block * t = (block *) malloc(sizeof(block));
        t->data = (double *) malloc ((n+1) * sizeof(double));
        t->size = n;
        return t;
     	}

void block_free(block * t) {
        free(t->data);
        free(t);
     	}

void block_copy(void *inp, void *outp) {
       	int i;
       	block * in = (block *) inp;
       	block * out = (block *) outp;

       	for(i=0; i< in->size; i++){
               out->data[i] = in->data[i];
       	}
       	out->size = in->size;
       	out->s = in-> s;
	}

void * block_copy_construct(void *xp) {
	block * x = (block *) xp;
	block * y = block_alloc(x->size);
	block_copy(x, y);
	return y;
        }

void block_destroy(void *xp){
	block_free( (block *) xp);
        }

///////////////////////////////Simulated annealing handling finished////////////

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

int deltaKronecker(const int i; const int j); 
int deltaKronecker(const int i; const int j) {
	int delta;
	if (i==j) delta = 1;
	else delta = 0;
	return delta;
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

//Simulated annealing functions
double L_function(void *xp)  
{
  block x = * ((block *) xp);
  //Data imports
  gsl_vector *saxs_ens = (gsl_vector *) (x->saxsEnsPtr);
  gsl_vector *saxs_exp = (gsl_vector *) (x->saxsExpPtr);
  gsl_vector *err_saxs = (gsl_vector *) (x->saxsErrPtr);
  gsl_matrix *saxs_pre = (gsl_matrix *) (x->saxsPrePtr);
  double saxs_scale = x->saxsScale;
  size_t L = x->size;
  size_t N = saxs_exp->size;

  double alpha_zero = 0.0;
  double Lfunc=0.0, fit_saxs=0.0, fit_saxs_mix = 0.0;
  for (int i = 0; i < L; i++)
	  alpha_zero+=x->alphas[i];

  Lfunc+=log(gsl_sf_gamma(alpha_zero)/gsl_sf_gamma(L/2));

  for (int i = 0; i < L; i++) {
	Lfunc+=log(gsl_sf_gamma(0.5)/gsl_sf_gamma(x->alphas[i]));
  }

  for (int i = 0; i < L; i++) {
        Lfunc+=(x->alphas[i]-0.5)*(gsl_sf_psi(x->alphas[i])-gsl_sf_psi(alpha_zero));
  }

  for( int i = 0; i< N; i++) { 
	fit_saxs += (pow( saxs_scale*gsl_vector_get(saxs_ens,i) - gsl_vector_get(saxs_exp,i), 2) / pow(gsl_vector_get(err_saxs,i),2) ); 
  }
  
  for( int i = 0; i< L; i++) {
	for (int j = 0; j < L; j++) {
		double mix_saxs = 0.0;	
		for (int k = 0; k < N; k++) {	 
        		mix_saxs += gsl_matrix_get(saxs_pre,k,i)*gsl_matrix_get(saxs_pre,k,j)/pow(gsl_vector_get(err_saxs,i),2);
		}
		fit_saxs_mix += mix_saxs * (x->alphas[i]*(alpha_zero - x->alphas[i])*deltaKronecker(i,j) + x->alphas[i]*x->alphas[j]*(1-deltaKronecker(i,j)));
	}
  } 
  fit_saxs_mix = 0.5 * fit_saxs_mix/(pow(alpha_zero,2)*(alpha_zero+1));
  L_func+=fit_saxs+fit_saxs_mix;
  return L_func;
}

double L_distance(void *xp, void *yp)
{
  block *x = ((block *) xp);
  block *y = ((block *) yp);
  double vector_distance = 0.0;
  for (int i=0; i<x->size; i++) {
	vector_distance+=gsl_pow_2(x->alphas[i]-y->alphas[i]);
  return sqrt(vector_distance);
}

void L_print (void *xp)
{
  block *x = (block *) xp;
  for(int i=0;i < x->size; i++){
      printf("%6.2lf ", x->data[i]);
  }
}

void L_take_step(const gsl_rng * r, void *xp, double step_size)
{
  block * x = ((block *) xp);
  //The index of which alpha should be modified
  int i = (int) round(gsl_rng_uniform(r)*x->size);
  x->alphas[i] = GSL_MAX(0, x->alphas[i]+gsl_ran_gaussian_ziggurat(r,step_size));
  //printf("%.2g\n", gsl_vector_get(x,i));

}


/*Overall algorithm
1. Read experimental data and parameter priors
2. Run simulated anealing to minimize function 
3. Iteratively remove structures with weights lower than wcut
*/
int main()
{
	//////////////////// Init section /////////////////////////////////////
	int n,k,steps,equilibration,samples;
	char mdfile[80], outfile[80];
	//Number of measurements in first curve
	int N; 
	char presaxsfile[80], saxsfile[80], saxserrfile[80]; 
	//Restart from already precaclculated vaules
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

	//Running params
	fscanf(stdin, "%s", &outfile); 
	fscanf(stdin, "%d", &equilibration); 
	fscanf(stdin, "%d", &steps); 
	fscanf(stdin, "%d", &samples); 

	double alpha_zero;

	gsl_matrix *saxs_pre = gsl_matrix_alloc(N,k);

	gsl_vector *saxs_exp = gsl_vector_alloc(N),
		*err_saxs = gsl_vector_alloc(N),
		*w_pre = gsl_vector_alloc(k),
		*w_ens_current = gsl_vector_alloc(k),
		*alpha_ens_current = gsl_vector_alloc(k),
		*saxs_ens_current = gsl_vector_alloc(N),
		*bayesian_weight1 = gsl_vector_alloc(k),
		*memory = gsl_vector_alloc(k+2);

	//Marks indexes that don't pass threshold filter
	bool removed_indexes[k];

	// Read in data from files //
	FILE * inFile = fopen(presaxsfile,"r"); gsl_matrix_fscanf(inFile,saxs_pre);fclose(inFile);
	inFile = fopen(saxsfile,"r"); gsl_vector_fscanf(inFile,saxs_exp); fclose(inFile);
	inFile = fopen(saxserrfile,"r"); gsl_vector_fscanf(inFile,err_saxs); fclose(inFile);
	inFile = fopen(mdfile,"r"); gsl_vector_fscanf(inFile,w_pre); fclose(inFile);

	// initialize random number generators //
	const gsl_rng_type *K; 
	gsl_rng *r; 
	gsl_rng_env_setup(); 
	K = gsl_rng_default;
	r = gsl_rng_alloc(K); 
	gsl_rng_set(r,time(NULL)); 

	block simAnBlock = block_alloc(k);

	//Initialize alphas with prior values
	for (int i = 0; i < N; i++) {
		simAnBlock->alphas[i] = gsl_vector_get(w_pre,i);
	}

	simAnBlock->saxsExpPtr = saxs_exp;
	simAnBlock->saxsErrPtr = err_saxs;
	simAnBlock->saxsEnsPtr = saxs_ens;
	simAnBlock->saxsPrePtr = saxs_pre;
	simAnBlock->saxsScale = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
	///////////////////////////////////////////////////////////////////////
	
	////////////////////// First iteration ////////////////////////////////
	cout<<"Equilibration started..."<<std::endl;
	int N_TRIES = 100*N;
	int ITERS_FIXED_T = 1;
	double STEP_SIZE = 1/N_TRIES;
	double K = 1/5;
	double T_INITIAL = 2; 
	double T_MIN = 2.0e-6;
	double MU_T = 1.003; 
	gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};
	//Define params before equilibration and after for next rounds
	gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, L_print,
		 block_copy, block_copy_construct, block_destroy,                
                 0, params);

	alpha_zero = 0.0;
	for (int i=0; i < k; i++) {
		alpha_zero+=simAnBlock->alphas[i];
		gsl_vector_set(alpha_ens_current,i,simAnBlock->alphas[i]);
	}
	for (int i=0; i < k; i++) {
                gsl_vector_set(w_ens_current,i,simAnBlock->alphas[i]/alpha_zero);
        }
	block_destroy(simAnBlock);
	/////////////////////////////////////////////////////////////////////
	
	///////////////////Next iterations //////////////////////////////////
	cout<<"Simulated annealing started"<<std::endl;
	int overall_iteration = 0;
	int L = k;
	int l;
	while ( gsl_vector_max(w_ens_current) < w_cut ) {
		
		block simAnBlock = block_alloc(L);
		gsl_matrix *saxs_pre_round = gsl_matrix_alloc(N,L);
		
		l = 0;
		for (int i = 0; i < k; i++) {
			if (removed_indexes[i]!=True) {
				for (int j = 0; j < N; j++) {
					gsl_matrix_set(saxs_pre_round,j,i);
				}
                		simAnBlock->alphas[l] = gsl_vector_get(alpha_ens_current,l);
				l++;
			}
        	}
		//saxs_exp and err_saxs are independent of run
        	simAnBlock->saxsExpPtr = saxs_exp;
        	simAnBlock->saxsErrPtr = err_saxs;
		simAnBlock->saxsPrePtr = saxs_pre_round;
		gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre_round, w_ens, 0.0, saxs_ens_current);
        	simAnBlock->saxsEnsPtr = saxs_ens_current;
		simAnBlock->saxsScale = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
		
		int N_TRIES = 100*L;
        	int ITERS_FIXED_T = 1; 
        	double STEP_SIZE = 1/N_TRIES;
        	double K = 1/5;
        	double T_INITIAL = 1;
        	double T_MIN = 2.0e-6;
        	double MU_T = 1.003;
        	gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

		//alphas are used from the previous simulation 
		gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, L_print,
                  block_copy, block_copy_construct, block_destroy, 
                  0, params);

		l=0;
		for (int i = 0; i < k; i++) if (removed_indexes[i] == True) l++;

		for (int i = 0; i < L; i++) alpha_zero +=simAnBlock->alphas[i];
		for (int i = 0; i < L; i++) {
			double wib = simAnBlock->alphas[i]/alpha_zero;
			if (wib < w_cut) {
				removed_indexes[i+l] = True;
			} else {
				gsl_vector_set(w_ens_current,i+l,wib);
			}
		}
	
		if (energyCurrent < energyMin) {
			energyMin = energyCurrent;
			last_updated = overall_iteration;
		}

		block_destroy(simAnBlock);	
		//If L_function doesn't improve after 10 iterations exit program
		if ((overall_iteration-last_updated)>10) break;
		overall_iteration++;
	}	
	///////////////////////////////////////////////////////////////////////	

	for( int l = 0; l < k; l++) {
        	gsl_vector_set(memory, l , gsl_vector_get(w_ens_current,l));
        }

	gsl_vector_set(memory, k, saxs_scale_current);
	gsl_vector_set(memory, k+1, energy_current);

	ofstream output(outfile);
        for (int i = 0; i < samples; i++) { for( int j = 0; j < n_sets*k+(2+n_sets); j++) { output << gsl_matrix_get(memory,i,j) << " "; if (j == n_sets*k+(1+n_sets)) { output << endl; } } }
        output.close();

	
	gsl_rng_free (r);
	return 0;
}
