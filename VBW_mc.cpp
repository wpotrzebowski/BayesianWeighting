#include "VBW_mc.hh"

using namespace std;
const double pi = M_PI;

block * block_alloc(size_t n) {
        block * t = (block *) malloc(sizeof(block));
        t->alphas = (double *) malloc ((n+1) * sizeof(double));
	//t->alphas = (double *) malloc ((n) * sizeof(double));
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
	//out->saxsEnsPtr = in->saxsEnsPtr;
	out->saxsPrePtr = in->saxsPrePtr;
	out->saxsMixPtr = in->saxsMixPtr;
       	out->OligomericSpecies = in->OligomericSpecies;
	out->Concentration = in->Concentration;
	out->MonomerMass = in->MonomerMass;
	out->OligomerOrder = in->OligomerOrder;
	out->saxsScale = in->saxsScale;	
	out->numberProcs = in->numberProcs;
	out->numberOfCurves = in->numberOfCurves;
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



void polySolver (int order, double cmass_ratio, gsl_vector *oligomeric_states, gsl_matrix *kconsts,  double *roots)
{

        /*Constans nomenclature
        k[0][0] - monomer always 1
        k[0][1] = [Dim]/[Mon]
        k[0]2] = [Tri]/[Mon]
        This will be expressed as monomer and different Ks will have to be recalculated
        */
        double *coefficents = (double * ) malloc( order * sizeof( double ));
        for (int i = 0; i<order; i++) {
                coefficents[i]=0;
        }

        //All fractions sum up to 1 (non x element in polynomial equation)
        coefficents[0]=-1;
        //Setup kconsts matrix
        for (int i = 1; i<order; i++) {
        	//if in equilibrium
        	if ( gsl_vector_get(oligomeric_states,i-1)>0.0 )
                	coefficents[i]=i*pow(cmass_ratio,i-1)*gsl_matrix_get(kconsts,0,i-1);
        }

        gsl_poly_complex_workspace * w
        = gsl_poly_complex_workspace_alloc( order );
        gsl_poly_complex_solve( coefficents, order, w, roots );
        gsl_poly_complex_workspace_free( w );
	free( coefficents );
}

/*void find_poly_root(gsl_vector *w_ens, gsl_vector *w_ens_prim, double ct, double ct_prim,
        double monomerMass, int k, int order, gsl_vector *oligomeric_species )
{
	double w_ens_[k];
	double w_ens_prim_[k];
	for (int i=0; i < k; i++) w_ens_[i] = gsl_vector_get(w_ens,i);
	find_poly_root(w_ens, w_ens_prim, ct, ct_prim,monomerMass, k, order, *oligomeric_species )
	for (int i=0; i < k; i++) w_ens_prim_[i] = gsl_vector_get(w_ens_prim,i);
}*/

void find_poly_root( gsl_vector *w_ens, gsl_vector *w_ens_prim, double ct, double ct_prim, 
	double monomerMass, int k, int order, gsl_vector *oligomeric_species )
{

	//This fucntion reads in sampled weights and initial concentration and
	// returns coupled weights given corresponding concentration
	//order is maximum oligomeric state order
	double *roots = (double * ) malloc( 2*(order-1) * sizeof( double ));
	//Monomers fraction in initial concentration
	//OBS!: We assume here that monomer will be first on the strucrture list
	double fm = gsl_vector_get(w_ens,0);
	double fm_prim;
	double cmass_ratio = monomerMass/ct;
	double cmass_ratio_prim = monomerMass/ct_prim;
	double cmass_ratio_prim_inv = ct_prim/monomerMass;
	double mono_fract;
	int N;
	double kN;
	double w;
	//Sum of K constants stored for given oligomeric specie
	gsl_matrix *KSums = gsl_matrix_alloc(order-2,order-1);
	//K consdtant stored for individual model
	gsl_vector *KConsts = gsl_vector_alloc(k-1);
	gsl_vector *oligomeric_states =  gsl_vector_alloc(order-1);

	for (int i = 0; i<order-1; i++) {
        	gsl_vector_set(oligomeric_states,i,0);
	}

	for (int i = 0; i<order-2; i++) {
                for (int j = 0; j<order-1; j++) {
                        gsl_matrix_set(KSums,i,j,0);
                }
        }
	//For monomer monomer is set to 1
	gsl_matrix_set(KSums,0,0,1);
	
	for(int i = 0; i < k; i++) {
		N = gsl_vector_get(oligomeric_species,i);
		if (N == 1) {
                        gsl_vector_set(oligomeric_states,0,1);
                }
                else {
			//Oligomeric state says one when there is given state
			gsl_vector_set(oligomeric_states,N-1,1);
			mono_fract = N*pow(fm,N); 
			w =  gsl_vector_get(w_ens,i);
			//TODO: Will have to something smarter
			//if (w < 0.001 ) w= 0.001;
                	kN = w*pow(cmass_ratio,N-1)/mono_fract;
			//kconsts are set with the resepect to monomer but it can be generalized
			gsl_vector_set(KConsts,i-1,kN);
                	gsl_matrix_set(KSums,0,N-1,gsl_matrix_get(KSums,0,N-1)+kN);
		}
        }


	polySolver(order,cmass_ratio_prim_inv,oligomeric_states,KSums,roots);
	
	//Output has to be reporocessed and wens_prum has to be updated
	//TODO: What if real non-negative solution is not found?
	for (int i = 0; i < order-1; i++)
    	{
      		if ((roots[2*i+1]) == 0.0 && roots[2*i]>0.0) {
			fm_prim = roots[2*i];
		}
	}
	gsl_vector_set(w_ens_prim, 0, fm_prim);
	for(int i = 1; i < k; i++) {
		N = gsl_vector_get(oligomeric_species,i);
                mono_fract = N*pow(fm_prim,N);
		gsl_vector_set(w_ens_prim,i,gsl_vector_get(KConsts,i-1)*mono_fract*pow(cmass_ratio_prim_inv,N-1));
	}
	free(roots);
	gsl_matrix_free(KSums);
	gsl_vector_free(KConsts);
	gsl_vector_free(oligomeric_states);	
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
  double monomerMass = x->MonomerMass;
  size_t nprocs = x->numberProcs;
  size_t Ncurves = x->numberOfCurves;
  size_t oligoOrder = x->OligomerOrder;
  size_t sample_size = x->size;
  double fit_saxs=0.0, fit_saxs_mix=0.0;

  //TODO: saxs_ens may not be necessary
  double *mix_saxs = (double *) ( x->saxsMixPtr );
  gsl_matrix *saxs_exp = (gsl_matrix *) ( x->saxsExpPtr );
  gsl_matrix *err_saxs = (gsl_matrix *) ( x->saxsErrPtr ); 
  //gsl_vector *saxs_ens = (gsl_vector *) (x->saxsEnsPtr);
  gsl_vector *saxs_scale = (gsl_vector *) (x->saxsScale);
  gsl_vector *oligomeric_species = (gsl_vector *) (x->OligomericSpecies);
  gsl_vector *concentration = (gsl_vector *) (x->Concentration);
  gsl_matrix *saxs_pre = (gsl_matrix *) (x->saxsPrePtr);
  size_t N = saxs_exp->size1;
  size_t L = saxs_pre->size2;
  //double *w_current = ( double * ) malloc (Ncurves * L * sizeof(double));
  //double *alpha_l = ( double * ) malloc ( N * sizeof(double));
  gsl_vector *w_current = gsl_vector_alloc(L);
  gsl_vector *w_current_prim = gsl_vector_alloc(L);
  gsl_vector *saxs_weights_ens = gsl_vector_alloc(N);
  gsl_vector *alpha_l = gsl_vector_alloc(L); 
  int rep = 0;
  double alpha_zero = 0.0;
  double log_gamma_2 = gsl_sf_lngamma(0.5);
  double Lfunc=0.0;
 
  for (int i = 0; i < L; i++) {
	  alpha_zero+=x->alphas[i];
  }

  Lfunc+= ( gsl_sf_lngamma(alpha_zero)-gsl_sf_lngamma(L/2) );

  for (int i = 0; i < L; i++) {
	Lfunc+=(log_gamma_2 - gsl_sf_lngamma( x->alphas[i] ));
  }

  for (int i = 0; i < L; i++) {
        Lfunc+=((x->alphas[i]-0.5)*(gsl_sf_psi(x->alphas[i])-gsl_sf_psi(alpha_zero)));
  	gsl_vector_set(w_current,i,x->alphas[i]/alpha_zero);
  }

  for( int l = 0; l < Ncurves; l++) {
	
	if (l == 0) { 
		//for( int i = 0; i< L; i++) alpha_l[i] = x->alphas[i];
		for( int i = 0; i< L; i++) gsl_vector_set(alpha_l,i, x->alphas[i]);
		gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_current, 0.0, saxs_weights_ens); 
	}
	else {
		find_poly_root(w_current, w_current_prim, gsl_vector_get(concentration,0), gsl_vector_get(concentration,l),
        	monomerMass, L, oligoOrder,oligomeric_species);
		//alpha_zero is sampled in each concentration 
		alpha_zero = x->alphas[L+l-1];
		//cout<<"Alpha zero for round "<<l<<" "<<alpha_zero<<std::endl;
		//for( int i = 0; i< L; i++) alpha_l[i] = gsl_vector_get(w_current_prim,i)*alpha_zero;
		for( int i = 0; i< L; i++) gsl_vector_set(alpha_l,i,gsl_vector_get(w_current_prim,i)*alpha_zero);
		gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_current_prim, 0.0, saxs_weights_ens);
	}

	for (int i = 0; i < N; i++) {
		fit_saxs += ( pow(gsl_vector_get(saxs_weights_ens,i) - gsl_matrix_get(saxs_exp,i,l), 2) / pow(gsl_matrix_get(err_saxs,i,l),2) );
	}
  	double smix, deltamix;
  	int i_ind,j_ind;
   
  	//gettimeofday(&t1, NULL);
  	/*#pragma omp parallel for \
  	default(none) shared(Ncurves, L, l, x, mix_saxs, alpha_l, alpha_zero, w_current)\
  	private (i_ind, j_ind, smix, deltamix) \
  	num_threads(nprocs) \
  	schedule(dynamic,16) \
  	reduction(+:fit_saxs_mix)*/

  	for(i_ind = 0; i_ind < L; i_ind++) {
  		for ( j_ind = i_ind; j_ind < L; j_ind++) {
			smix = mix_saxs[ Ncurves*L*l + L*i_ind + j_ind ];
			//deltamix = (i_ind!=j_ind) ? -2*alpha_l[i_ind]*alpha_l[j_ind] : alpha_l[i_ind]*(alpha_zero - alpha_l[i_ind]);
			deltamix = (i_ind!=j_ind) ? -2*gsl_vector_get(alpha_l,i_ind)*gsl_vector_get(alpha_l,j_ind) :\
				 gsl_vector_get(alpha_l,i_ind)*(alpha_zero - gsl_vector_get(alpha_l,i_ind));
               		fit_saxs_mix += deltamix * smix;
       		}
  	}
  	//gettimeofday(&t2, NULL);
  	//TODO: Alpha zero parameters will be sampled extra  
  
  	fit_saxs_mix /= (pow(alpha_zero,2)*(alpha_zero+1));
  	Lfunc+=0.5*(fit_saxs+fit_saxs_mix);
  }
  // compute and print the elapsed time in millisec
  //elapsedTime = (t2.tv_sec - t1.tv_sec)*1000.0;      // sec to ms
  //elapsedTime += (t2.tv_usec - t1.tv_usec)/1000.0;
  //cout << "Time: "<< fit_saxs_mix<< " : "<<elapsedTime << " ms."<<std::endl;
  gsl_vector_free(saxs_weights_ens);
  gsl_vector_free(w_current);
  gsl_vector_free(w_current_prim);
  gsl_vector_free(alpha_l);

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
  size_t Ncurves = x->numberOfCurves;
  gsl_matrix *saxs_pre = (gsl_matrix *) (x->saxsPrePtr);
  size_t L = saxs_pre->size2;
  for(int i=0; i < L; i++){
  	alpha_zero += x->alphas[i];
  }
  cout<<" Weights: ";
  for(int i=0; i < L; i++){
      weight =  x->alphas[i]/alpha_zero;
      cout<<weight<<" ";
     //Add vector save here
  }
  cout<<" Alphas zeros: ";
  for(int l = 0; l<Ncurves - 1; l++) 
	cout<<x->alphas[L+l]<<" ";
  cout<<std::endl;
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
	const int &N, const std::string &presaxsfile, const int &Ncurves, const std::string &curvesfile,
	const std::string &filelist, const double &monomerMass,
	const std::string &outfile, const int &nprocs, const double &w_cut)
{
	//////////////////// Init section /////////////////////////////////////
	//Number of saxs measurements is read from individual curves
	//double wdelta = 0.0001;	
	int read_success =0;
	gsl_siman_params_t params;
	int N_TRIES; //Seems to be inactive?
        int ITERS_FIXED_T ;
        double STEP_SIZE;
        double K;
        double T_INITIAL;
        double MU_T;
        double T_MIN;

	/*Input Format of main file
	Not stdin but read main file
	1. Start form saved state? 0 - no, 1 - yes
	2. File containing list of structural models in format "oligomeric_state,filename", e.g. "2,dimer_1.pdb"
	3. File with prior weights. Order is the same as for list of structures 
	4. File contatiming theorrtcial scattering intesities [matrix]
	5. File with scattering curves in format "qvector intensity error"
	6. Output file
	7. Number of processors used in simulations
	8. Weights cuts
	*/

	double alpha_zero;
	double energy_current, energy_min;
	float acceptance_rate = 1.0;
	size_t Ncurvesm1;
	
	double *saxs_mix;
        saxs_mix =  (double * ) malloc( Ncurves * k * k * sizeof( double ));
	if (saxs_mix==NULL) { 
		cerr<<"Cannot allocate memory"<<std::endl;
		exit(1);
	}
	gsl_vector *saxs_scale_current = gsl_vector_alloc(Ncurves), 
		*alpha_ens_current = gsl_vector_alloc(k),
		*tostart = gsl_vector_alloc(k+2),
		*memory = gsl_vector_alloc(k+2),
		*bayesian_weight1 = gsl_vector_alloc(k),
                *bayesian_weight1_current = gsl_vector_alloc(k);

	//Done to satisfy non-zero vector lenght
	if (Ncurves>1)  Ncurvesm1 = Ncurves-1; 
	else  Ncurvesm1 = 1;  
	gsl_vector *alpha_zeros_prim =  gsl_vector_alloc(Ncurvesm1);
 
	gsl_matrix *saxs_exp = gsl_matrix_alloc(N,Ncurves),
                *err_saxs = gsl_matrix_alloc(N,Ncurves);
	
	gsl_vector *saxs_exp_vec = gsl_vector_alloc(N),
                *err_saxs_vec = gsl_vector_alloc(N);

	gsl_vector *w_current[Ncurves],
                *w_ens_current[Ncurves],
		*saxs_ens_current[Ncurves];

	gsl_vector *oligomeric_species = gsl_vector_alloc(k);
        gsl_vector *concentrations = gsl_vector_alloc(Ncurves);
        int oligomerOrder=0;//Maximum value from oligomeric_species + 1 (due to polysolver definition);

	gsl_matrix *saxs_pre = gsl_matrix_alloc(N,k);
	//SAXS file matrix containing q vectors, Intenisty, erros
	gsl_matrix *saxs_file_matrix = gsl_matrix_alloc(N,3);

	for (int i = 0; i < Ncurves; i++) {
		w_current[i] = gsl_vector_alloc(k);
		w_ens_current[i] = gsl_vector_alloc(k);	
		saxs_ens_current[i] = gsl_vector_alloc(N);
	}
	gsl_vector_set_zero(bayesian_weight1);
	//TODO: Samples, set to maximum 500, which is also the maximum number of iterations.
	int samples = 500;
	gsl_matrix *weight_samples = gsl_matrix_alloc(samples,k);;

	//Marks indexes that don't pass threshold filter
	bool removed_indexes[k];
	for (int i = 0; i < k; i++) removed_indexes[i]=false;

	// Prior file is read //
	FILE* inFile = fopen(mdfile.c_str(),"r"); gsl_vector_fscanf(inFile,w_current[0]); fclose(inFile);
	cout<<"Loaded priors"<<std::endl;	
	//Reading file list
	std::ifstream inFileList(filelist.c_str());
        int list_line = 0;
	int oligomeric_state;
	std::string structure_file;
	while (inFileList >> oligomeric_state >> structure_file) {
		if (oligomeric_state > oligomerOrder) 
			oligomerOrder = oligomeric_state;
		gsl_vector_set(oligomeric_species,list_line,oligomeric_state);
		list_line++;
	}	
	//+1 beacuse of poly solver convention
	oligomerOrder +=1;	

	if (list_line != k) {
                cerr<<"Number of records in file list doesn't agree with sumber of simualted curves"<<std::endl;
                exit (EXIT_FAILURE);
        }

	cout<<"Loaded file list"<<std::endl;	
	//Experimental SAXS files are read
	std::ifstream inSAXSfile(curvesfile.c_str());
	int curves_file_line = 0;
	float conc; 
	std::string saxs_dat_file;
	while (inSAXSfile >> conc >> saxs_dat_file) {
		gsl_vector_set(concentrations,curves_file_line,conc);
		FILE *inSAXSdat = fopen(saxs_dat_file.c_str(),"r");
		gsl_matrix_fscanf(inSAXSdat,saxs_file_matrix); 
		for (int i = 0;  i< N; i++) {
			gsl_matrix_set(saxs_exp,i,curves_file_line,gsl_matrix_get(saxs_file_matrix,i,1));
			gsl_matrix_set(err_saxs,i,curves_file_line,gsl_matrix_get(saxs_file_matrix,i,2));
		} 
		fclose(inSAXSdat);
		curves_file_line++;
	}
	
	//fclose(inSAXSfile);
	if (curves_file_line != Ncurves) {
		cerr<<"Incompatible number of scattering curves"<<std::endl;
		exit (EXIT_FAILURE);
	}	

	cout<<"Loaded experimental file"<<std::endl;
	//Simulated SAXS data is read
	//TODO: The assumption is that the number of points is the same all experimental ana simulated	
	inFile = fopen(presaxsfile.c_str(),"r"); gsl_matrix_fscanf(inFile,saxs_pre);fclose(inFile);
	cout<<"Files reading finished"<<std::endl;

	//TODO: This will be in NCurves dimenssion

	// initialize random number generators //
	const gsl_rng_type *Krng; 
	gsl_rng *r; 
	gsl_rng_env_setup(); 
	Krng = gsl_rng_default;
	r = gsl_rng_alloc(Krng); 
	gsl_rng_set(r,time(NULL)); 
	
	//Ncurves -1 accomodates extra parameters needed for coupling
	block *simAnBlock = block_alloc(k + Ncurves - 1);
	
	//Initialize alphas with prior values
	for (int i = 0; i < k; i++) {
		simAnBlock->alphas[i] = gsl_vector_get(w_current[0],i);
	}
	//Initializing alphas_zero with 1.0
	for (int l = 0; l < Ncurves-1; l++) {
                simAnBlock->alphas[k+l] = 1.0;
        }
	simAnBlock->saxsExpPtr = saxs_exp;
	simAnBlock->saxsErrPtr = err_saxs;
	simAnBlock->saxsPrePtr = saxs_pre;
	simAnBlock->numberProcs = nprocs;
	simAnBlock->numberOfCurves = Ncurves;
	for (int i=0; i<Ncurves; i++) {
		*saxs_exp_vec = gsl_matrix_column(saxs_exp,i).vector;
		*err_saxs_vec = gsl_matrix_column(err_saxs,i).vector;
		gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_current[i], 0.0, saxs_ens_current[i]);
		gsl_vector_set(saxs_scale_current, i, SaxsScaleMean(saxs_ens_current[i],\
			saxs_exp_vec, err_saxs_vec ,N));
	}
	simAnBlock->saxsScale = saxs_scale_current;
	//simAnBlock->saxsEnsPtr = saxs_ens_current;

	simAnBlock->OligomericSpecies = oligomeric_species;
        simAnBlock->Concentration = concentrations;
        simAnBlock->MonomerMass = monomerMass;
        simAnBlock->OligomerOrder = oligomerOrder;
		
	if(again == 1){ inFile = fopen("restart.dat","r"); gsl_vector_fscanf(inFile,tostart); fclose(inFile); }	
	//timeval t1, t2;
	//double elapsedTime;
  	//gettimeofday(&t1, NULL);

	//double *saxs_mix =  (double * ) malloc( Ncurves * k * k * sizeof( double ));
	double smix;
	//double cs_mix;
        //#pragma omp parallel for reduction(+:smix) reduction(+:cs_mix) num_threads(nprocs) 	
	for (int l = 0; l < Ncurves; l++) {
		//#pragma omp parallel for reduction(+:smix) num_threads(nprocs)    
		for( int i = 0; i< k; i++) {
        		for (int j = 0; j < k; j++) {
				smix = 0.0;
                		for (int m = 0; m < N; m++) {
					smix+=gsl_matrix_get(saxs_pre,m,i)*gsl_matrix_get(saxs_pre,m,j)/pow(gsl_matrix_get(err_saxs,m,l),2);   
				}
                        	saxs_mix[ Ncurves*k*l + k*i + j ] = smix;
        		}
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
		gsl_vector_set(saxs_scale_current, 0, gsl_vector_get(tostart,k+1));
		simAnBlock->saxsScale = saxs_scale_current;
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
	
		for (int j=0; j < Ncurves-1; j++) {
                        gsl_vector_set(alpha_zeros_prim,j,simAnBlock->alphas[k+j]);
			cout<<"Alphas zeros prim "<< gsl_vector_get(alpha_zeros_prim,j)<<std::endl;;
                }

		for (int i=0; i < k; i++) {
                	gsl_vector_set(w_ens_current[0],i,gsl_vector_get(alpha_ens_current,i)/alpha_zero);
        	}
		
		energy_min = L_function(simAnBlock);
		for (int i= 0; i < Ncurves; i++) {
			if ( i > 0 ) { 
				//TODO: Do it proper - repetition
				find_poly_root(w_ens_current[0], w_ens_current[i], gsl_vector_get(concentrations,0), gsl_vector_get(concentrations,i),
                		monomerMass, k, oligomerOrder,oligomeric_species);
			}
			*saxs_exp_vec = gsl_matrix_column(saxs_exp,i).vector;
	                *err_saxs_vec = gsl_matrix_column(err_saxs,i).vector;
        	        gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_ens_current[i], 0.0, saxs_ens_current[i]);
              		gsl_vector_set(saxs_scale_current, i, SaxsScaleMean(saxs_ens_current[i],\
                        saxs_exp_vec, err_saxs_vec ,N));
		}
		block_destroy(simAnBlock);
		free(saxs_mix);
		/////////////////////////////////////////////////////////////////////
		//Store alphas after equilibration stage
		ofstream restart("restart.dat");
        	for(int j = 0; j < k; j++) { restart << gsl_vector_get(alpha_ens_current,j)<<" "; }
        	restart <<energy_min<<" "<<gsl_vector_get(saxs_scale_current,0)<<std::endl;
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
		
		gsl_matrix *saxs_pre_round = gsl_matrix_alloc(N,L);
                double  *saxs_mix_round =  (double * ) malloc( Ncurves * L * L * sizeof( double ));
		block *simAnBlock = block_alloc(L + Ncurves - 1 );
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
		for (int j = 0; j < Ncurves-1; j++) {
			simAnBlock->alphas[L+j] = gsl_vector_get(alpha_zeros_prim,j);
		}
		for (int c=0; c < Ncurves; c++) {
			//#pragma omp parallel for reduction(+:smix) num_threads(nprocs) 
                	for( int i = 0; i < L; i++) {
                        	for (int j = 0; j < L; j++) {
					smix = 0.0;
                                	for (int m = 0; m < N; m++) {
						smix+=gsl_matrix_get(saxs_pre_round,m,i)*gsl_matrix_get(saxs_pre_round,m,j)/pow(gsl_matrix_get(err_saxs,m,c),2);
        	                	}
					saxs_mix_round[ Ncurves*L*c + L*i + j ] = smix;
                       		}
			}
                }
		cout<<"Molecular mass"<<monomerMass<<" "<<oligomerOrder<<std::endl;
		//saxs_exp and err_saxs are independent of run
        	simAnBlock->saxsExpPtr = saxs_exp;
        	simAnBlock->saxsErrPtr = err_saxs;
		simAnBlock->saxsPrePtr = saxs_pre_round;
        	simAnBlock->saxsMixPtr = saxs_mix_round;
		//simAnBlock->saxsEnsPtr = saxs_ens_current;
		simAnBlock->saxsScale = saxs_scale_current;
		simAnBlock->numberProcs = nprocs;	
		simAnBlock->numberOfCurves = Ncurves;
	
		simAnBlock->OligomericSpecies = oligomeric_species;
	        simAnBlock->Concentration = concentrations;
        	simAnBlock->MonomerMass = monomerMass;
        	simAnBlock->OligomerOrder = oligomerOrder;	
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

		cout<<"SimAn finsihed "<<oligomerOrder<<std::endl;
		energy_current = L_function(simAnBlock);

                //If L_function doesn't improve after 10 iterations exit program
		newL = 0;
		m = 0; 
		alpha_zero = 0.0;
		for ( int i = 0; i < L; i++ ) alpha_zero +=simAnBlock->alphas[i];

		//Sampled alpha zeros are stored here
		for ( int j = 0; j < Ncurves-1; j++ ) gsl_vector_set(alpha_zeros_prim,j,simAnBlock->alphas[L+j]);
		
		double new_alpha_zero = 0.0;
		for ( int i = 0; i < k; i++ ) {
			if ( removed_indexes[i]==false ) {
				double wib = simAnBlock->alphas[m]/alpha_zero;
				if ( wib < w_cut ) {
					gsl_vector_set( alpha_ens_current, i, 0.0 );
					gsl_vector_set( w_ens_current[0], i, 0.0);
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
				gsl_vector_set( w_ens_current[0],i,gsl_vector_get(alpha_ens_current,i)/new_alpha_zero );
			}
		}
		//Stoping simulations if weights don't change for more than delta (0.001)
		//if (wdelta_count == newL) {cout<<"Simulations stopped because weights don't progress"<<std::endl; break;}
		for ( int i = 0; i < Ncurves; i++ ) {	
			for (int j = 0; j < L; j++) cout<<gsl_vector_get(oligomeric_species,j);
			cout<<" Concentrations: "<<gsl_vector_get(concentrations,i)<<std::endl;
			if ( i > 0 ) {
                                //TODO: Do it proper. Also oligomric speciec will have to change
                                find_poly_root(w_ens_current[0], w_ens_current[i], gsl_vector_get(concentrations,0), gsl_vector_get(concentrations,i),
                                monomerMass, L, oligomerOrder,oligomeric_species);
                        }
			for (int j = 0; j < L; j++) cout<<gsl_vector_get(w_ens_current[0],j)<<" "<<gsl_vector_get(w_ens_current[i],j)<<std::endl;
			*saxs_exp_vec = gsl_matrix_column(saxs_exp,i).vector;
	                *err_saxs_vec = gsl_matrix_column(err_saxs,i).vector;
        	        gsl_blas_dgemv(CblasNoTrans, 1.0, saxs_pre, w_ens_current[i], 0.0, saxs_ens_current[i]);
                	gsl_vector_set(saxs_scale_current, i, SaxsScaleMean(saxs_ens_current[i],\
                        	saxs_exp_vec, err_saxs_vec ,N));
		}
		//gsl_blas_dgemv(CblasNoTrans, 1.0, cs_pre, w_ens_current, 0.0, cs_ens_current);	
		//Structural library size after discarding structures with weight lower than cuttof
		L = newL;
		cout<<"Before destroy"<<std::endl;
		block_destroy(simAnBlock);	
		overall_iteration++;
		cout<<"After destroy"<<std::endl;	
		if (energy_current < energy_min) {	
			energy_min = energy_current;
                        last_updated = overall_iteration;

			for( int l = 0; l < k; l++) {
                		gsl_vector_set(memory, l , gsl_vector_get(w_ens_current[0],l));
        		}

        		gsl_vector_set(memory, k, gsl_vector_get(saxs_scale_current,0));
        		gsl_vector_set(memory, k+1, energy_current);

			ofstream output(outfile,  std::ofstream::out | std::ofstream::trunc);
        		//All weights plus saxs scale factor
        		for( int j = 0; j < k + 1; j++) output << gsl_vector_get(memory,j) << " ";
        		output <<gsl_vector_get(memory,k+1)<<endl;
			output.close();
		}
	
		sampling_step = overall_iteration-1;
		for (int jind=0; jind<k; jind++) {
                    gsl_matrix_set(weight_samples,sampling_step,jind,gsl_vector_get(w_ens_current[0],jind));
                }

                double niter = 1.0/double(sampling_step+1);
               	gsl_vector_add(bayesian_weight1,w_ens_current[0]);
                gsl_vector_memcpy(bayesian_weight1_current,bayesian_weight1);
               	gsl_vector_scale(bayesian_weight1_current,niter);

		free(saxs_mix_round);
                gsl_matrix_free(saxs_pre_round);

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
	for (int i = 0; i < Ncurves; i++) {
		gsl_vector_free(saxs_ens_current[i]);
		gsl_vector_free(w_ens_current[i]);
		gsl_vector_free(w_current[i]);
	}
	free(saxs_ens_current);
	free(w_ens_current);
	free(w_ens_current);
	free(saxs_mix);

}

