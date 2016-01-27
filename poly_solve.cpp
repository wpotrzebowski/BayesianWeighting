#include <stdio.h>
#include <cmath>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <iostream>

using namespace std;
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

        for(int i = 1; i < k; i++) {
                N = gsl_vector_get(oligomeric_species,i);
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

	for (int i = 0; i<order-1; i++) {
		cout<<gsl_matrix_get(KSums,0,i)<<std::endl;
	}
        polySolver(order,cmass_ratio_prim_inv,oligomeric_states,KSums,roots);

        //Output has to be reporocessed and wens_prum has to be updated
        //TODO: What if real non-negative solution is not found?
        for (int i = 0; i < order-1; i++)
        {
                if ((roots[2*i+1]) == 0.0 && roots[2*i]>0.0 ) {
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

int main(void)
{


double ct = 0.25;
double ct_prim = 0.5;
//double ctot2 = 0.5;
double monomerMass = 24.14E+3;
int k = 3;
int order = 5;

gsl_vector *w_ens =  gsl_vector_alloc(3);
gsl_vector *w_ens_prim =  gsl_vector_alloc(3);
gsl_vector *oligomeric_species =  gsl_vector_alloc(3);

gsl_vector_set(w_ens,0,0.87);
gsl_vector_set(w_ens,1,0.12);
gsl_vector_set(w_ens,2,0.01);


gsl_vector_set(oligomeric_species,0,1);
gsl_vector_set(oligomeric_species,1,2);
gsl_vector_set(oligomeric_species,2,4);

find_poly_root( w_ens, w_ens_prim, ct, ct_prim, monomerMass, k, order, oligomeric_species );
for (int i = 0; i<3; i++) 
	cout<<gsl_vector_get(w_ens,i)<<" "<<gsl_vector_get(w_ens_prim,i)<<std::endl;

return 0;
}

