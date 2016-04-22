#include "VBW_CS_sc.hh"
#include "LFuncGpu.hh"

using namespace std;
const double pi = M_PI;

block * block_alloc(size_t n) {
        block * t = (block *) malloc(sizeof(block));
        t->alphas = (double *) malloc ((n+1) * sizeof(double));
        t->size = n;
        return t;
     	}

void block_free(block * t) {
        free(t->alphas);
        free(t);
     	}

void block_copy(void *inp, void *outp) {
       	int i;
       	block * in = (block *) inp;
       	block * out = (block *) outp;

       	for(i=0; i< in->size; i++){
               out->alphas[i] = in->alphas[i];
       	}
       	out->size = in->size;
	out->saxsExpPtr = in->saxsExpPtr;
	out->saxsErrPtr = in->saxsErrPtr;
	//TODO: Seems not to be used
	out->saxsEnsPtr = in->saxsEnsPtr;
	out->saxsPrePtr = in->saxsPrePtr;
	out->saxsMixPtr = in->saxsMixPtr;
       	out->saxsScale = in->saxsScale;
	
	out->csExpPtr = in->csExpPtr;
        out->csErrPtr = in->csErrPtr;
	out->csRmsPtr = in->csRmsPtr;
	//TODO: Seems not to be used
        out->csEnsPtr = in->csEnsPtr;
        out->csPrePtr = in->csPrePtr;
        out->csMixPtr = in->csMixPtr;
	
	out->numberProcs = in->numberProcs;
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

double ientropy(const gsl_vector *w, int i) {
	double ie = 0.0;
	//for (int i=0; i<k; i++) 
	ie = gsl_vector_get(w,i)*log2(gsl_vector_get(w,i));
	return ie;
}

double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k) {

	double jsd=0.0, s1=0.0, s2=0.0;
	for (int i=0; i<k; i++) {
		if ( gsl_vector_get(w_a,i) == 0.0 || gsl_vector_get(w_b,i) == 0.0) continue;
		s1 +=  gsl_vector_get(w_a,i)*log2(2*gsl_vector_get(w_a,i)/(gsl_vector_get(w_a,i)+gsl_vector_get(w_b,i)));
		s2 +=  gsl_vector_get(w_b,i)*log2(2*gsl_vector_get(w_b,i)/(gsl_vector_get(w_a,i)+gsl_vector_get(w_b,i)));
	}
	jsd =  0.5*(s1+s2);
	return jsd;
}

//TODO: These will be replaced with other function
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

double SaxsScaleMean(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N)
{
	double tempa = 0.0, tempb = 0.0;
	for( int i = 0; i< N; i++) {
		tempa += gsl_vector_get(saxs_ens,i)*gsl_vector_get(saxs_exp,i)/gsl_vector_get(err_saxs,i);
		tempb += pow(gsl_vector_get(saxs_ens,i),2.0)/gsl_vector_get(err_saxs,i);
	}
	return tempa/tempb;
}

double SaxsScaleStandardDeviation(gsl_vector *saxs_ens, gsl_vector *saxs_exp, gsl_vector *err_saxs, int N, double T)
{
	double temp = 0.0;
	for( int i = 0; i< N; i++) { 
		temp += pow(gsl_vector_get(saxs_ens,i),2.0)/gsl_vector_get(err_saxs,i); 
	}
	return sqrt(T/temp);
}

///////////////////Simulated annealing functions////////////////////////////////
double L_function(void *xp)  
{
  //timeval t1, t2;
  //double elapsedTime;
  //gettimeofday(&t1, NULL);

  block *x = (block *) xp;
  //Data imports
  //gsl_vector *saxs_ens = (gsl_vector *) (x->saxsEnsPtr);
  gsl_vector *saxs_exp = (gsl_vector *) (x->saxsExpPtr);
  gsl_vector *err_saxs = (gsl_vector *) (x->saxsErrPtr);
  gsl_matrix *saxs_pre = (gsl_matrix *) (x->saxsPrePtr);
  
  //gsl_vector *cs_ens = (gsl_vector *) (x->csEnsPtr);
  gsl_vector *cs_exp = (gsl_vector *) (x->csExpPtr);
  gsl_vector *cs_err = (gsl_vector *) (x->csErrPtr);
  gsl_vector *cs_rms = (gsl_vector *) (x->csRmsPtr);
  gsl_matrix *cs_pre = (gsl_matrix *) (x->csPrePtr);

  double *mix_saxs = (double *) (x->saxsMixPtr);
  double *mix_cs = (double *) (x->csMixPtr);
  double saxs_scale = x->saxsScale;
  int nprocs = x->numberProcs;
  size_t L = x->size;
  size_t N = saxs_exp->size;
  size_t n = cs_exp->size;
  int rep = 0;
  double alpha_zero = 0.0;
  double saxs_alpha_ens[N];
  double cs_alpha_ens[n];
  double log_gamma_2 = gsl_sf_lngamma(0.5);
  double Lfunc=0.0;
  double fit_saxs=0.0, fit_saxs_mix = 0.0;
  double fit_cs=0.0, fit_cs_mix = 0.0;
  
  for (int i = 0; i < L; i++)
	  alpha_zero+=x->alphas[i];

  Lfunc+= ( gsl_sf_lngamma(alpha_zero)-gsl_sf_lngamma(L/2) );

  for (int i = 0; i < L; i++) {
	Lfunc+=(log_gamma_2 - gsl_sf_lngamma( x->alphas[i] ));
  }

  for (int i = 0; i < L; i++) {
        Lfunc+=((x->alphas[i]-0.5)*(gsl_sf_psi(x->alphas[i])-gsl_sf_psi(alpha_zero)));
  }


  for( int i = 0; i< N; i++) {
	saxs_alpha_ens[i] = 0.0;
	for (int k = 0; k < L; k++) {
		saxs_alpha_ens[i]+=gsl_matrix_get(saxs_pre,i,k)*x->alphas[k];
	}
	fit_saxs += ( pow(saxs_alpha_ens[i]/alpha_zero - gsl_vector_get(saxs_exp,i), 2) / pow(gsl_vector_get(err_saxs,i),2) );
  }
 
  for( int i = 0; i< n; i++) {
        cs_alpha_ens[i] = 0.0;
        for (int k = 0; k < L; k++) {
                cs_alpha_ens[i]+=gsl_matrix_get(cs_pre,i,k)*x->alphas[k];
        }
        fit_cs += ( pow(cs_alpha_ens[i]/alpha_zero - gsl_vector_get(cs_exp,i), 2) / ( pow(gsl_vector_get(cs_err,i),2) + pow(gsl_vector_get(cs_rms,i),2) ) );
  }



  double smix, csmix, deltamix;
  int i_ind,j_ind;
   
  //gettimeofday(&t1, NULL);
  #pragma omp parallel for \
  default(none) shared(L,x,mix_saxs,mix_cs,alpha_zero,nprocs)\
  private (i_ind, j_ind, smix, csmix, deltamix) \
  num_threads(nprocs) \
  schedule(dynamic,nprocs) \
  reduction(+:fit_saxs_mix)\
  reduction(+:fit_cs_mix)

  for(i_ind = 0; i_ind < L; i_ind++) {
  	for ( j_ind = i_ind; j_ind < L; j_ind++) {
		smix = mix_saxs[L*i_ind+j_ind];
		csmix = mix_cs[L*i_ind+j_ind];
                deltamix = (i_ind!=j_ind) ? -2*x->alphas[i_ind]*x->alphas[j_ind] : x->alphas[i_ind]*(alpha_zero - x->alphas[i_ind]);
               	fit_saxs_mix += deltamix * smix;
		fit_cs_mix += deltamix * csmix;
       }
  }
  //gettimeofday(&t2, NULL); 
  fit_saxs_mix /= (pow(alpha_zero,2)*(alpha_zero+1));
  fit_cs_mix /= (pow(alpha_zero,2)*(alpha_zero+1));
  //Lfunc+=0.5*(fit_saxs+fit_saxs_mix);
  Lfunc+=0.5*(fit_saxs+fit_cs+fit_saxs_mix+fit_cs_mix);

  // compute and print the elapsed time in millisec
  //elapsedTime = (t2.tv_sec - t1.tv_sec)*1000.0;      // sec to ms
  //elapsedTime += (t2.tv_usec - t1.tv_usec)/1000.0;
  //cout << "Time: "<< fit_saxs_mix<< " : "<<elapsedTime << " ms."<<std::endl;

  return Lfunc;
}

double L_distance(void *xp, void *yp)
{
  block *x = (block *) xp;
  block *y = (block *) yp;
  double vector_distance = 0.0;
  for (int i=0; i<x->size; i++) {
	vector_distance+=gsl_pow_2(x->alphas[i]-y->alphas[i]);
  }
  return sqrt(vector_distance);
}

//No printing is done by default
void L_print (void *xp)
{
  block *x = (block *) xp;
  double alpha_zero = 0.0;
  double weight; 
  for(int i=0; i < x->size; i++){
  	alpha_zero += x->alphas[i];
  }
  for(int i=0; i < x->size; i++){
      weight =  x->alphas[i]/alpha_zero;
     //Add vector save here
  }
}

void L_take_step(const gsl_rng * r, void *xp, double step_size)
{
  	block * x = (block *) xp;
  	//The index of which alpha should be modified
  	int i = (int) round(gsl_rng_uniform(r)*x->size);
  	double u = x->alphas[i]+gsl_ran_gaussian_ziggurat(r,step_size);
	x->alphas[i] = GSL_MAX(0.001, u); 
}


/*Overall algorithm
1. Read experimental data and parameter priors
2. Run simulated anealing to minimize function 
3. Iteratively remove structures with weights lower than wcut
*/
void run_vbw(const int &again, const int &k, const std::string &mdfile,
        const int &N, const int &n, const int &Ncurves,
        const std::string &presaxsfile, const std::string &saxsfile,
        const std::string &precsfile, const std::string &csfile,
        const std::string &outfile, const int &nprocs, const double &w_cut)
{
	//////////////////// Init section /////////////////////////////////////
	double saxs_scale_current;
	double wdelta = 0.0001;

	gsl_siman_params_t params;
	int N_TRIES; //Seems to be inactive?
    	int ITERS_FIXED_T ;
    	double STEP_SIZE;
    	double K;
    	double T_INITIAL;
    	double MU_T;
    	double T_MIN;


	double alpha_zero;
	double energy_current, energy_min;
	double *saxs_mix; 
	double *cs_mix;
	float acceptance_rate = 1.0;
 	saxs_mix = (double * ) malloc( k * k * sizeof( double ));
	cs_mix = (double * ) malloc( k * k * sizeof( double ));
	gsl_matrix *saxs_pre = gsl_matrix_alloc(N,k);
	gsl_matrix *cs_pre = gsl_matrix_alloc(n,k);
   	gsl_matrix *saxs_file_matrix = gsl_matrix_alloc(N,3);
    	gsl_matrix *cs_file_matrix = gsl_matrix_alloc(n,3);

	gsl_vector *saxs_exp = gsl_vector_alloc(N),
		*err_saxs = gsl_vector_alloc(N),
		*cs_exp = gsl_vector_alloc(n),
        	*cs_err = gsl_vector_alloc(n),
		*cs_rms = gsl_vector_alloc(n),
		*w_pre = gsl_vector_alloc(k),
		*w_ens_current = gsl_vector_alloc(k),
		*alpha_ens_current = gsl_vector_alloc(k),
		*tostart = gsl_vector_alloc(k+2),
		*saxs_ens_current = gsl_vector_alloc(N),
		*cs_ens_current = gsl_vector_alloc(n),
		*memory = gsl_vector_alloc(k+2),
		*bayesian_weight1 = gsl_vector_alloc(k),
        	*bayesian_weight1_current = gsl_vector_alloc(k);

	gsl_vector_set_zero(bayesian_weight1);
	//TODO: Samples, set to maximum 500, which is also the maximum number of iterations.
	int samples = 500;
	gsl_matrix *weight_samples = gsl_matrix_alloc(samples,k);;

	//Marks indexes that don't pass threshold filter
	bool removed_indexes[k];
	for (int i = 0; i < k; i++) removed_indexes[i]=false;

    //Read prior files
	FILE *inFile = fopen(mdfile.c_str(),"r");
	gsl_vector_fscanf(inFile,w_pre); fclose(inFile);

	// Read in data from files //
    inFile = fopen(presaxsfile.c_str(),"r");
    gsl_matrix_fscanf(inFile,saxs_pre); fclose(inFile);

    inFile = fopen(precsfile.c_str(),"r");
    gsl_matrix_fscanf(inFile,cs_pre); fclose(inFile);

	//Read scattering file
    FILE *inSAXSdat = fopen(saxsfile.c_str(),"r");
    gsl_matrix_fscanf(inSAXSdat,saxs_file_matrix);
    for (int i = 0;  i< N; i++) {
        gsl_vector_set(saxs_exp,i,gsl_matrix_get(saxs_file_matrix,i,1));
       	gsl_vector_set(err_saxs,i,gsl_matrix_get(saxs_file_matrix,i,2));
    }
    fclose(inSAXSdat);

	FILE *inCSdat = fopen(csfile.c_str(),"r");
    gsl_matrix_fscanf(inCSdat,cs_file_matrix);
    for (int i = 0;  i< N; i++) {
        gsl_vector_set(cs_exp,i,gsl_matrix_get(cs_file_matrix,i,1));
       	gsl_vector_set(cs_err,i,gsl_matrix_get(cs_file_matrix,i,2));
    }
    fclose(inCSdat);

	cout<<"Files reading finished"<<std::endl;

	// initialize random number generators //
	const gsl_rng_type *Krng; 
	gsl_rng *r; 
	gsl_rng_env_setup(); 
	Krng = gsl_rng_default;
	r = gsl_rng_alloc(Krng); 
	gsl_rng_set(r,time(NULL)); 

	block *simAnBlock = block_alloc(k);
	
	//Initialize alphas with prior values
	for (int i = 0; i < k; i++) {
		simAnBlock->alphas[i] = gsl_vector_get(w_pre,i);
	}
	simAnBlock->saxsExpPtr = saxs_exp;
	simAnBlock->saxsErrPtr = err_saxs;
	simAnBlock->saxsPrePtr = saxs_pre;
	simAnBlock->csExpPtr = cs_exp;
    	simAnBlock->csErrPtr = cs_err;
	simAnBlock->csRmsPtr = cs_rms;
    	simAnBlock->csPrePtr = cs_pre;
	simAnBlock->numberProcs = nprocs;
	gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_pre, 0.0, saxs_ens_current);	
	gsl_blas_dgemv(CblasNoTrans, 1.0, cs_pre, w_pre, 0.0, cs_ens_current);
	saxs_scale_current = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
	simAnBlock->saxsScale = saxs_scale_current;
	simAnBlock->saxsEnsPtr = saxs_ens_current;
	simAnBlock->csEnsPtr = cs_ens_current;
	
	if(again == 1){ inFile = fopen("restart.dat","r"); gsl_vector_fscanf(inFile,tostart); fclose(inFile); }	
	//timeval t1, t2;
	//double elapsedTime;
  	//gettimeofday(&t1, NULL);

	double smix;
	double csmix;
    	#pragma omp parallel for reduction(+:smix) reduction(+:csmix) num_threads(nprocs)
	//#pragma omp parallel for reduction(+:smix) num_threads(nprocs)    
	for( int i = 0; i< k; i++) {
        	for (int j = 0; j < k; j++) {
			smix = 0.0;
			csmix = 0.0;
                	for (int m = 0; m < N; m++) {
				smix+=gsl_matrix_get(saxs_pre,m,i)*gsl_matrix_get(saxs_pre,m,j)/pow(gsl_vector_get(err_saxs,m),2);   
			}
			for (int m = 0; m < n; m++) {
                                csmix+=gsl_matrix_get(cs_pre,m,i)*gsl_matrix_get(cs_pre,m,j)/(pow(gsl_vector_get(cs_err,m),2)+pow(gsl_vector_get(cs_rms,m),2));
                        }
                        saxs_mix[i*k+j] = smix;
			cs_mix[i*k+j] = csmix;
        	}
	}
	/*gettimeofday(&t2, NULL);
	// compute and print the elapsed time in millisec
	elapsedTime = (t2.tv_sec - t1.tv_sec)*1000.0;      // sec to ms
  	elapsedTime += (t2.tv_usec - t1.tv_usec)/1000.0;
  	cout << "Time: "<< elapsedTime << " ms."<<std::endl;*/

	simAnBlock->saxsMixPtr = saxs_mix;
	simAnBlock->csMixPtr = cs_mix;
	///////////////////////////////////////////////////////////////////////

	cout<<"Values have been set"<<std::endl;
	///////////////////////////////////////////////////////////////////////
	if(again == 1) { 
		inFile = fopen("restart.dat","r"); 
		gsl_vector_fscanf(inFile,tostart); 
		fclose(inFile); 
		for( int i = 0; i< k; i++) gsl_vector_set(alpha_ens_current,i,gsl_vector_get(tostart,i));
		energy_min = gsl_vector_get(tostart,k);
		simAnBlock->saxsScale = gsl_vector_get(tostart,k+1);
	}
	else {	
		////////////////////// First iteration ////////////////////////////////
		cout<<"Equilibration started..."<<std::endl;
		
		N_TRIES = 1; //Seems to be inactive?
		ITERS_FIXED_T = 1;
		STEP_SIZE = 1;
		K = 1.0;
		T_INITIAL = 2.0; 
		MU_T = 1.000025;
       		T_MIN = 2.7776e-11;
		params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

		//Define params before equilibration and after for next rounds
		gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, NULL,
		 	block_copy, block_copy_construct, block_destroy,                
                 	0, params, &acceptance_rate);

		alpha_zero = 0.0;
		for (int i=0; i < k; i++) {
			alpha_zero+=simAnBlock->alphas[i];
			gsl_vector_set(alpha_ens_current,i,simAnBlock->alphas[i]);
		}
		for (int i=0; i < k; i++) {
                	gsl_vector_set(w_ens_current,i,gsl_vector_get(alpha_ens_current,i)/alpha_zero);
        	}
		energy_min = L_function(simAnBlock);
		gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_ens_current, 0.0, saxs_ens_current);
		saxs_scale_current = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
		gsl_blas_dgemv(CblasNoTrans, 1.0, cs_pre, w_ens_current, 0.0, cs_ens_current);
		block_destroy(simAnBlock);
		free(saxs_mix);
		free(cs_mix);
		/////////////////////////////////////////////////////////////////////
	
		//Store alphas after equilibration stage
		ofstream restart("restart.dat");
        	for(int j = 0; j < k; j++) { restart << gsl_vector_get(alpha_ens_current,j)<<" "; }
        	restart <<energy_min<<" "<<saxs_scale_current<<std::endl;
        	restart.close();
	}
			
	///////////////////Next iterations //////////////////////////////////
	cout<<"Simulated annealing started"<<std::endl;
	int overall_iteration = 0;
	int sampling_step;
	int last_updated;
	int L = k;
	int l, m, newL;
	//Energy from first iteration
	while ( L > 1 ) {
		cout<<"Starting "<<overall_iteration+1<<" iteration with "<<L<<" models"<<std::endl;
		
		block *simAnBlock = block_alloc(L);
		gsl_matrix *saxs_pre_round = gsl_matrix_alloc(N,L);
		double  *saxs_mix_round =  (double * ) malloc( k * k * sizeof( double )); 

		gsl_matrix *cs_pre_round = gsl_matrix_alloc(n,L);
                double  *cs_mix_round =  (double * ) malloc( k * k * sizeof( double ));

		l = 0;
		for (int i = 0; i < k; i++) {
			if (removed_indexes[i]==false) {
				for (int j = 0; j < N; j++) {
					gsl_matrix_set(saxs_pre_round,j,l,gsl_matrix_get(saxs_pre,j,i));
				}
				for (int j = 0; j < n; j++) {
                                        gsl_matrix_set(cs_pre_round,j,l,gsl_matrix_get(cs_pre,j,i));
                                }
                		simAnBlock->alphas[l] = gsl_vector_get(alpha_ens_current,i);
				l++;
			}
        	}
		//#pragma omp parallel for reduction(+:smix) reduction(+:cs_mix) num_threads(nprocs)  
		#pragma omp parallel for reduction(+:smix) reduction(+:csmix) num_threads(nprocs) 
                for( int i = 0; i < L; i++) {
                        for (int j = 0; j < L; j++) {
				smix = 0.0;
				csmix = 0.0;
                                for (int m = 0; m < N; m++) {
					smix+=gsl_matrix_get(saxs_pre_round,m,i)*gsl_matrix_get(saxs_pre_round,m,j)/pow(gsl_vector_get(err_saxs,m),2);
        	                }
				for (int m = 0; m < n; m++) {
                                        csmix+=gsl_matrix_get(cs_pre_round,m,i)*gsl_matrix_get(cs_pre_round,m,j)/(pow(gsl_vector_get(cs_err,m),2)+pow(gsl_vector_get(cs_rms,m),2));
                                }
	                        saxs_mix_round[i*L+j]=smix;
				cs_mix_round[i*L+j]=csmix;
                       	}
		}
                
		//saxs_exp and err_saxs are independent of run
        	simAnBlock->saxsExpPtr = saxs_exp;
        	simAnBlock->saxsErrPtr = err_saxs;
		simAnBlock->saxsPrePtr = saxs_pre_round;
        	simAnBlock->saxsMixPtr = saxs_mix_round;
		simAnBlock->saxsEnsPtr = saxs_ens_current;
		simAnBlock->saxsScale = saxs_scale_current;

		simAnBlock->csExpPtr = cs_exp;
                simAnBlock->csErrPtr = cs_err;
		simAnBlock->csRmsPtr = cs_rms;
                simAnBlock->csPrePtr = cs_pre_round;
                simAnBlock->csMixPtr = cs_mix_round;
                simAnBlock->csEnsPtr = cs_ens_current;

		simAnBlock->numberProcs = nprocs;	
		
		////////////////////////Short equilibration period to find step size/////////////////////////
		N_TRIES = 1;
                ITERS_FIXED_T = 1000;
                K = 1.0;
                T_INITIAL = 1.0;
                MU_T = 1.00005;
                T_MIN = 1.0;
		//Itertate over different step size
		float dmin = 10;
		for (double s=0.01; s<2.1; s+=0.1) { 
                	params = {N_TRIES, ITERS_FIXED_T, s, K, T_INITIAL, MU_T, T_MIN};
		
                	//alphas are used from the previous simulation 
                	gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, NULL,
                                block_copy, block_copy_construct, block_destroy,
                                0, params, &acceptance_rate);
			if(fabs(acceptance_rate -0.5) < dmin) { 
				dmin = fabs(acceptance_rate -0.5);
				STEP_SIZE = s;		
			}	
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		cout<<"STEP_SIZE set to: "<<STEP_SIZE<<std::endl;
		N_TRIES = 1;
        	ITERS_FIXED_T = 1; 
        	STEP_SIZE = 1;
        	K = 1.0;
        	T_INITIAL = 1.0;
		MU_T = 1.00005;
		T_MIN = 1.3888e-11;
        	params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

		//alphas are used from the previous simulation 
		gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, NULL,
                  		block_copy, block_copy_construct, block_destroy, 
                  		0, params, &acceptance_rate);

		energy_current = L_function(simAnBlock);

                //If L_function doesn't improve after 10 iterations exit program
		newL = 0;
		m = 0; 
		alpha_zero = 0.0;
		for ( int i = 0; i < L; i++ ) alpha_zero +=simAnBlock->alphas[i];

		double new_alpha_zero = 0.0;
		for ( int i = 0; i < k; i++ ) {
			if ( removed_indexes[i]==false ) {
				double wib = simAnBlock->alphas[m]/alpha_zero;
				if ( wib < w_cut ) {
					gsl_vector_set( alpha_ens_current, i, 0.0 );
					gsl_vector_set( w_ens_current, i, 0.0);
					removed_indexes[i] = true;
				} else {
					new_alpha_zero += simAnBlock->alphas[m];
					gsl_vector_set( alpha_ens_current, i, simAnBlock->alphas[m] );
					newL++;
				}
				m++;
			}
		}
		
		//int wdelta_count = 0;
		for ( int i = 0; i < k; i++ ) {
                        if (removed_indexes[i]==false) {
				gsl_vector_set( w_ens_current,i,gsl_vector_get(alpha_ens_current,i)/new_alpha_zero );
			}
		}
		//Stoping simulations if weights don't change for more than delta (0.001)
		//if (wdelta_count == newL) {cout<<"Simulations stopped because weights don't progress"<<std::endl; break;}

		gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_ens_current, 0.0, saxs_ens_current);
	        saxs_scale_current = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
		gsl_blas_dgemv(CblasNoTrans, 1.0, cs_pre, w_ens_current, 0.0, cs_ens_current);	
		//Structural library size after discarding structures with weight lower than cuttof
		L = newL;
		
		block_destroy(simAnBlock);	
		overall_iteration++;
	
		if (energy_current < energy_min) {	
			energy_min = energy_current;
                        last_updated = overall_iteration;

			for( int l = 0; l < k; l++) {
                		gsl_vector_set(memory, l , gsl_vector_get(w_ens_current,l));
        		}

        		gsl_vector_set(memory, k, saxs_scale_current);
        		gsl_vector_set(memory, k+1, energy_current);

			ofstream output(outfile,  std::ofstream::out | std::ofstream::trunc);
        		//All weights plus saxs scale factor
        		for( int j = 0; j < k + 1; j++) output << gsl_vector_get(memory,j) << " ";
        		output <<gsl_vector_get(memory,k+1)<<endl;
			output.close();
		}
	
		sampling_step = overall_iteration-1;
		for (int jind=0; jind<k; jind++) {
                    gsl_matrix_set(weight_samples,sampling_step,jind,gsl_vector_get(w_ens_current,jind));
                }

                double niter = 1.0/double(sampling_step+1);
               	gsl_vector_add(bayesian_weight1,w_ens_current);
                gsl_vector_memcpy(bayesian_weight1_current,bayesian_weight1);
               	gsl_vector_scale(bayesian_weight1_current,niter);

		free(saxs_mix_round);
		gsl_matrix_free(saxs_pre_round);
		free(cs_mix_round);
                gsl_matrix_free(cs_pre_round);
		if ((overall_iteration-last_updated)>10) {
                        cout<<"Energy hasn't decreased for 10 iterations. Stopping simulations"<<std::endl;
                        break;
                }
		
		if (overall_iteration == samples) {
			cout<<"Maximum number of iteration has been reached. Stopping simulation"<<std::endl;
			break;
		}

	}	
	///////////////////////////////////////////////////////////////////////	

	//Calculating posterior expected divergence
        //TODO: Make a cluean-up with vector
        double jsd1_sum = 0.0;
        double jsd1 = 0.0;
        for (int s=0; s<sampling_step; s++) {
                for (int j=0; j<k; j++) {
                        gsl_vector_set(bayesian_weight1,j,gsl_matrix_get(weight_samples,s,j));
                }
                jsd1 = jensen_shannon_div(bayesian_weight1_current,bayesian_weight1,k);
                jsd1_sum += sqrt(jsd1);
        }
        cout<<"\nPED1: "<<jsd1_sum/double(sampling_step)<<" from "<<sampling_step<<" steps"<<std::endl;

	gsl_rng_free (r);
}
