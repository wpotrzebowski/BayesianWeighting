#include "VBW_sc.hh"
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
	out->saxsEnsPtr = in->saxsEnsPtr;
	out->saxsPrePtr = in->saxsPrePtr;
	out->saxsMixPtr = in->saxsMixPtr;
       	out->saxsScale = in->saxsScale;
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

double ientropy(const gsl_vector *w, int k) {
	double ie = 0.0;
	for (int i=0; i<k; i++) 
		ie -= gsl_vector_get(w,i)*log2(gsl_vector_get(w,i));
	return ie;
}

double jensen_shannon_div(const gsl_vector *w_a, const gsl_vector *w_b, int k) {

	double jsd;
	gsl_vector *w_c = gsl_vector_alloc(k);
	gsl_vector_memcpy(w_c, w_b);
	gsl_vector_add(w_c,w_a);
	gsl_vector_scale(w_c,0.5);
	jsd = ientropy(w_c,k) - 0.5*ientropy(w_a,k) - 0.5*ientropy(w_b,k);
	return jsd;
}


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
  /*timeval t1, t2;
  double elapsedTime;
  gettimeofday(&t1, NULL);*/

  block *x = (block *) xp;
  //Data imports
  gsl_vector *saxs_ens = (gsl_vector *) (x->saxsEnsPtr);
  gsl_vector *saxs_exp = (gsl_vector *) (x->saxsExpPtr);
  gsl_vector *err_saxs = (gsl_vector *) (x->saxsErrPtr);
  gsl_matrix *saxs_pre = (gsl_matrix *) (x->saxsPrePtr);
  double *mix_saxs = (double *) (x->saxsMixPtr);
  double saxs_scale = x->saxsScale;
  int nprocs = x->numberProcs;
  size_t L = x->size;
  size_t N = saxs_exp->size;
  int rep = 0;
  double alpha_zero = 0.0;
  double alpha_ens[N];
  double log_gamma_2 = gsl_sf_lngamma(0.5);
  double Lfunc=0.0, fit_saxs=0.0, fit_saxs_mix = 0.0;

  
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
	alpha_ens[i] = 0.0;
	for (int k = 0; k < L; k++) {
		alpha_ens[i]+=gsl_matrix_get(saxs_pre,i,k)*x->alphas[k];
	}
	fit_saxs += ( pow(alpha_ens[i]/alpha_zero - gsl_vector_get(saxs_exp,i), 2) / pow(gsl_vector_get(err_saxs,i),2) );
  }

  //gettimeofday(&t1, NULL);

  double smix, deltamix;
  int i_ind,j_ind;
  #pragma omp parallel for reduction(+:fit_saxs_mix) num_threads(nprocs)
  for( i_ind = 0; i_ind < L; i_ind++) {
  	for (j_ind = i_ind; j_ind < L; j_ind++) {
        	smix = mix_saxs[L*i_ind+j_ind];
		//if (i_ind!=j_ind) {
		//      	deltamix = - x->alphas[i_ind]*x->alphas[j_ind];
                //} else {
		// 	deltamix = x->alphas[i_ind]*(alpha_zero - x->alphas[i_ind]);
                //}
                deltamix = (i_ind!=j_ind) ? -2*x->alphas[i_ind]*x->alphas[j_ind] : x->alphas[i_ind]*(alpha_zero - x->alphas[i_ind]);
               	fit_saxs_mix += smix * deltamix;
       }
  }
  //gettimeofday(&t2, NULL); 
  fit_saxs_mix /= (pow(alpha_zero,2)*(alpha_zero+1));
  Lfunc+=0.5*(fit_saxs+fit_saxs_mix);
  // compute and print the elapsed time in millisec
  /*elapsedTime = (t2.tv_sec - t1.tv_sec)*1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec)/1000.0;
  cout << "Time: "<< elapsedTime << " ms."<<std::endl;*/

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

void L_print (void *xp)
{
  block *x = (block *) xp;
  double alpha_zero = 0.0;
  int weights_over_t=0;
  double xdata; 
  for(int i=0; i < x->size; i++){
  	alpha_zero += x->alphas[i];
  }
  for(int i=0; i < x->size; i++){
      xdata =  x->alphas[i]/alpha_zero;
      if (xdata > 0.01) weights_over_t++;
      //printf("%6.2lf ", xdata);
  }
  printf("Weights over threshold = %d ", weights_over_t);
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
int main()
{
	//////////////////// Init section /////////////////////////////////////
	int n,k,nprocs,again=0;
	double saxs_scale_current;
	float w_cut;
	char mdfile[80], outfile[80];
	int N; 
	char presaxsfile[80], saxsfile[80], saxserrfile[80]; 
	double wdelta = 0.0001;	
	int read_success =0;

	read_success = fscanf(stdin, "%d", &again); //Should precompuated values be used?
	read_success = fscanf(stdin, "%d", &k); //Number of structures 
	read_success = fscanf(stdin, "%s", &mdfile[0]); //Prior weights
	read_success = fscanf(stdin, "%d", &N); //Number of SAXS measurments
	read_success = fscanf(stdin, "%s", &presaxsfile[0]); //Theoretical SAXS curves
	read_success = fscanf(stdin, "%s", &saxsfile[0]); //Experimental SAXS files
	read_success = fscanf(stdin, "%s", &saxserrfile[0]); //Experimental Errors

	read_success = fscanf(stdin, "%s", &outfile[0]); 
	read_success = fscanf(stdin, "%d", &nprocs); 
	read_success = fscanf(stdin, "%f", &w_cut); //Weight cutoff

	if (read_success == 0) { 
		cerr<<"Error reading files"<<std::endl;
		exit (EXIT_FAILURE);
	}
	double alpha_zero;
	double energy_current, energy_min;
	double *saxs_mix; 
 	saxs_mix = (double * ) malloc( k * k * sizeof( double ));
	gsl_matrix *saxs_pre = gsl_matrix_alloc(N,k);

	gsl_vector *saxs_exp = gsl_vector_alloc(N),
		*err_saxs = gsl_vector_alloc(N),
		*w_pre = gsl_vector_alloc(k),
		*w_ens_current = gsl_vector_alloc(k),
		*alpha_ens_current = gsl_vector_alloc(k),
		*tostart = gsl_vector_alloc(k+1),
		*saxs_ens_current = gsl_vector_alloc(N),
		*memory = gsl_vector_alloc(k+2);

	//Marks indexes that don't pass threshold filter
	bool removed_indexes[k];
	for (int i = 0; i < k; i++) removed_indexes[i]=false;

	// Read in data from files //
	FILE * inFile = fopen(presaxsfile,"r"); gsl_matrix_fscanf(inFile,saxs_pre);fclose(inFile);
	inFile = fopen(saxsfile,"r"); gsl_vector_fscanf(inFile,saxs_exp); fclose(inFile);
	inFile = fopen(saxserrfile,"r"); gsl_vector_fscanf(inFile,err_saxs); fclose(inFile);
	inFile = fopen(mdfile,"r"); gsl_vector_fscanf(inFile,w_pre); fclose(inFile);
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
	simAnBlock->numberProcs = nprocs;
	gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_pre, 0.0, saxs_ens_current);	
	saxs_scale_current = SaxsScaleMean(saxs_ens_current,saxs_exp,err_saxs,N);
	simAnBlock->saxsScale = saxs_scale_current;
	simAnBlock->saxsEnsPtr = saxs_ens_current;
	
	if(again == 1){ inFile = fopen("restart.dat","r"); gsl_vector_fscanf(inFile,tostart); fclose(inFile); }	
	//timeval t1, t2;
	//double elapsedTime;
  	//gettimeofday(&t1, NULL);

	double smix;
        #pragma omp parallel for reduction(+:smix) num_threads(nprocs) 	
	for( int i = 0; i< k; i++) {
        	for (int j = 0; j < k; j++) {
			smix = 0.0;
                	for (int m = 0; m < N; m++) {
				smix+=gsl_matrix_get(saxs_pre,m,i)*gsl_matrix_get(saxs_pre,m,j)/pow(gsl_vector_get(err_saxs,m),2);   
			}
                        saxs_mix[i*k+j] = smix;
        	}
	}
	/*gettimeofday(&t2, NULL);
	// compute and print the elapsed time in millisec
	elapsedTime = (t2.tv_sec - t1.tv_sec)*1000.0;      // sec to ms
  	elapsedTime += (t2.tv_usec - t1.tv_usec)/1000.0;
  	cout << "Time: "<< elapsedTime << " ms."<<std::endl;*/

	simAnBlock->saxsMixPtr = saxs_mix;

	///////////////////////////////////////////////////////////////////////

	cout<<"Values have been set"<<std::endl;
	///////////////////////////////////////////////////////////////////////
	if(again == 1) { 
		inFile = fopen("restart.dat","r"); 
		gsl_vector_fscanf(inFile,tostart); 
		fclose(inFile); 
		for( int i = 0; i< k; i++) gsl_vector_set(alpha_ens_current,i,gsl_vector_get(tostart,i));
		energy_min = gsl_vector_get(tostart,k);
	}
	else {	
		////////////////////// First iteration ////////////////////////////////
		cout<<"Equilibration started..."<<std::endl;
		int N_TRIES = 1;
		int ITERS_FIXED_T = 1;
		double STEP_SIZE = 1;
		double K = 1.0;
		double T_INITIAL = 2.0; 
        	//double MU_T = 1.0105;
		double MU_T = 1.000025;
       		double T_MIN = 2.7776e-11;
		gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

		//Define params before equilibration and after for next rounds
		gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, NULL,
		 	block_copy, block_copy_construct, block_destroy,                
                 	0, params);

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
		block_destroy(simAnBlock);
		free(saxs_mix);
		/////////////////////////////////////////////////////////////////////
	
		//Store alphas after equilibration stage
		ofstream restart("restart.dat");
        	for(int j = 0; j < k; j++) { restart << gsl_vector_get(alpha_ens_current,j)<<" "; }
        	restart <<energy_min<<std::endl;
        	restart.close();
	}
			
	///////////////////Next iterations //////////////////////////////////
	cout<<"Simulated annealing started"<<std::endl;
	int overall_iteration = 0;
	int last_updated;
	int L = k;
	int l, m, newL;
	//Energy from first iteration
	while ( L > 0 ) {
		cout<<"Starting "<<overall_iteration+1<<" iteration with "<<L<<" models"<<std::endl;
		
		block *simAnBlock = block_alloc(L);
		gsl_matrix *saxs_pre_round = gsl_matrix_alloc(N,L);
		double  *saxs_mix_round =  (double * ) malloc( k * k * sizeof( double )); 
		l = 0;
		for (int i = 0; i < k; i++) {
			if (removed_indexes[i]==false) {
				for (int j = 0; j < N; j++) {
					gsl_matrix_set(saxs_pre_round,j,l,gsl_matrix_get(saxs_pre,j,i));
				}
                		simAnBlock->alphas[l] = gsl_vector_get(alpha_ens_current,i);
				l++;
			}
        	}

		#pragma omp parallel for reduction(+:smix) num_threads(nprocs) 
                for( int i = 0; i < L; i++) {
                        for (int j = 0; j < L; j++) {
				smix = 0.0;
                                for (int m = 0; m < N; m++) {
					smix+=gsl_matrix_get(saxs_pre_round,m,i)*gsl_matrix_get(saxs_pre_round,m,j)/pow(gsl_vector_get(err_saxs,m),2);
        	                }
	                        saxs_mix_round[i*L+j]=smix;
                       	}
		}
                
		//saxs_exp and err_saxs are independent of run
        	simAnBlock->saxsExpPtr = saxs_exp;
        	simAnBlock->saxsErrPtr = err_saxs;
		simAnBlock->saxsPrePtr = saxs_pre_round;
        	simAnBlock->saxsMixPtr = saxs_mix_round;
		simAnBlock->saxsEnsPtr = saxs_ens_current;
		simAnBlock->saxsScale = saxs_scale_current;
		simAnBlock->numberProcs = nprocs;	
	
		int N_TRIES = 1;
        	int ITERS_FIXED_T = 1; 
        	double STEP_SIZE = 1;
        	double K = 1.0;
        	double T_INITIAL = 1.0;
        	//double MU_T = 1.0211;
		double MU_T = 1.00005;
		double T_MIN = 1.3888e-11;
        	gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};

		//alphas are used from the previous simulation 
		gsl_siman_solve(r, simAnBlock, L_function, L_take_step, L_distance, NULL,
                  		block_copy, block_copy_construct, block_destroy, 
                  		0, params);
		/*cout<<"Current weights: ";
		for ( int i = 0; i < k; i++ ) {
			cout<< gsl_vector_get( w_ens_current,i )<<" ";
		}*/
		energy_current = L_function(simAnBlock);
                if (energy_current < energy_min) {
                        energy_min = energy_current;
                        last_updated = overall_iteration;
                } 
		else {
			cout<<"Energy increases... Stoping simulations"<<std::endl;
			break;
		}
                //If L_function doesn't improve after 10 iterations exit program
                if ((overall_iteration-last_updated)>10) break;

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
		
		//Structural library size after discarding structures with weight lower than cuttof
		L = newL;
		
		block_destroy(simAnBlock);	
		overall_iteration++;
		
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
		
		free(saxs_mix_round);
		gsl_matrix_free(saxs_pre_round);
	}	
	///////////////////////////////////////////////////////////////////////	

	gsl_rng_free (r);
	return 0;
}
