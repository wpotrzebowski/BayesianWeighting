#include <stdio.h>
#include <cmath>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void polySolver (int order, double cmass_ratio, gsl_vector *oligomeric_species, gsl_matrix *kconsts,  double roots[])
{

	/*Constans nomenclature
	k[0][0] - monomer always 1
	k[0][1] = [Dim]/[Mon]
	k[0]2] = [Tri]/[Mon]
	This will be expressed as monomer and different Ks will have to be recalculated
	*/
	double coefficents[order];
	for (int i = 0; i<order; i++) {
        	coefficents[i]=0;
	}


	//All fractions sum up to 1 (non x element in polynomial equation)
	coefficents[0]=-1;
	//Setup kconsts matrix
	for (int i = 1; i<order; i++) {
	//if in equilibrium
	if ( gsl_vector_get(oligomeric_species,i-1)>0.0 )
		coefficents[i]=i*pow(cmass_ratio,i-1)*gsl_matrix_get(kconsts,0,i-1);
	}

  	gsl_poly_complex_workspace * w
      	= gsl_poly_complex_workspace_alloc (order);

  	gsl_poly_complex_solve (coefficents, order, w, roots);

  	gsl_poly_complex_workspace_free (w);
}

int main(void)
{

int order=5;
double roots[2*(order-1)];
double coefficents[order];

double ctot1 = 0.25;
//double ctot2 = 0.5;
double molecularMass = 27.98E+3;
double cmass_ratio = ctot1/molecularMass;
double kdsum = 8.26E+3;
double ktsum = kdsum*kdsum*2.83E+2;

gsl_matrix *kconsts = gsl_matrix_alloc(order-2,order-1);
gsl_vector *oligomeric_species =  gsl_vector_alloc(order-1);

for (int i = 0; i<order-1; i++) {
        gsl_vector_set(oligomeric_species,i,0);
}

gsl_vector_set(oligomeric_species,0,1);
gsl_vector_set(oligomeric_species,1,1);
gsl_vector_set(oligomeric_species,3,1);

for (int i = 0; i<order-2; i++) {
        for (int j = 0; j<order-1; j++) {
                gsl_matrix_set(kconsts,i,j,1);
        }
}

gsl_matrix_set(kconsts,0,1,kdsum);
gsl_matrix_set(kconsts,0,3,ktsum);

polySolver(order,cmass_ratio,oligomeric_species,kconsts,roots);
printf ("Conc = %+.8f \n",ctot1);
for (int i = 0; i < order-1; i++)
    {
      if ((roots[2*i+1]) == 0.0 && roots[2*i]>0.0) {
      	printf (" %.18f",
              roots[2*i]);
      	printf (", %.18f",
              2*kdsum*pow(roots[2*i],2)*cmass_ratio);
      	printf (", %.18f\n",
              4*ktsum*pow(roots[2*i],4)*pow(cmass_ratio,3));
      }
    }

return 0;
}
