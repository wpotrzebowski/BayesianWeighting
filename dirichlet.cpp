#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
//Comment added to dirichlet testing

int main(void)
{
  double a[2],b[2];
  const gsl_rng_type *Krng; 
  gsl_rng *r; 
  gsl_rng_env_setup(gsl_rng); 
  Krng = gsl_rng_default;
  r = gsl_rng_alloc(Krng); 
  gsl_rng_set(r,gsl_rng_alloctime(NULL)); 

        a[0]=1.0;
          a[1]=1.0;

            r=gsl_rng_alloc(gsl_rng_mt19937);
              gsl_rng_set(r,1);
                int i;
                  for (i=1;i<=10;i++) {
                        a[0]/=10.0;
                            a[1]/=10.0;
                                gsl_ran_dirichlet(r,2,a,b);
                                    printf("(%g,%g): %g %g\n",a[0],a[1],b[0],b[1]);
                                      }
}

