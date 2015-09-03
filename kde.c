/**
 * @file kde.c
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
 * Estimates the kernel density p(x) at a given value x from
 * an array of sample points.  Uses the default algorithm from
 * the R langauge's 'density' function.  Requires the GSL statistics
 * library.  
 *   
 * @section LICENCE
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */
 
 
#include "kde.h"
 
/** Estimate bandwidth using Silverman's "rule of thumb" 
 * (Silverman 1986, pg 48 eq 3.31).  This is the default
 * bandwith estimator for the R 'density' function.  */
double nrd0(gsl_matrix *samples, const int J, const int N) {
	double * x;
	x = new double [N];
	for (int i=0; i<N; i++) x[i] = gsl_matrix_get(samples,i,J); 

	gsl_sort(x, 1, N);
	double hi = gsl_stats_sd(x, 1, N);
	double iqr =
		gsl_stats_quantile_from_sorted_data (x,1, N,0.75) - 
        gsl_stats_quantile_from_sorted_data (x,1, N,0.25);
	double lo = GSL_MIN(hi, iqr/1.34);
	double bw = 0.9 * lo * pow(N,-0.2);
	delete x;
	return(bw);
}
 
/* kernels for kernel density estimates */
double gauss_kernel(double x) { 
	return exp(-(gsl_pow_2(x)/2))/(M_SQRT2*sqrt(M_PI)); 
}
 
long double kerneldensity(gsl_matrix *samples, gsl_vector *obs, int n, int d) {
	int i, j;
	gsl_vector *h_vec = gsl_vector_alloc(d);
	long double sum = 0;
	double h; 
	long double prod;
	for(j=0; j < d; j++) { 
		h = GSL_MAX(nrd0(samples, j, n), 1e-6);
		gsl_vector_set(h_vec,j,h);
	}
	for(i=0; i < n; i++) {
		prod = 1.0;
		for(j=0; j < d; j++) {
			prod*=gauss_kernel((gsl_vector_get(obs,j) - gsl_matrix_get(samples,i,j))/gsl_vector_get(h_vec,j))/gsl_vector_get(h_vec,j);
		}
		sum += prod;
	}
	return sum/n;
}

