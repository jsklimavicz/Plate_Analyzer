/*
stats/funclib/cloglik.c
This file is a component of Plate_Analyzer, which can count Drosophila
L1 larvae, classify them as alive or dead, and determine dose-repsonse
information based on live/dead count data. 

Copyright (C) 2021 James Klimavicz

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "multimin.h"
#include "ll_mod.h"
#include "cloglik.h"
#include "ls.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <time.h>

void ll3_min(double *b, //Vector of parameters for optimization
	double *probs, //Array of probs
	double *conc, //Array of concentrations
	int probs_size, //Length of the prob/conc arrays
	double minimum, //The function value
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	double *beta, //2-value array with parameters for beta prior on b2
	const int method) //ID of method to be usedfor optimization. 
{
	struct ll_param fparam = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1  b2
	unsigned type[3] =  {3,    3,  3};

	multimin(3, b, &minimum, type, lb, ub, &ll3f,&ll3df,&ll3fdf,
			(void *) &fparam, optim_par);
}

void ll2_min(double *b, //Vector of parameters for optimization
	double *probs, //Array of probs
	double *conc, //Array of concentrations
	int probs_size, //Length of the prob/conc arrays
	double minimum, //The function value
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	const int method) //ID of method to be usedfor optimization.
{
	struct ll_param fparam = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};
	
	//                   b0	  b1 
	unsigned type[2] =  {3,    3};

	multimin(2, b, &minimum, type, lb, ub, &ll2f,&ll2df,&ll2fdf,
			(void *) &fparam, optim_par);
}

void ll3_array_min(
	int probs_size, //Number of vals in the prob/conc arrays
	int n_iters, //Number of different sets of probs
	double b[][3], //The vals to be optimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, //Array of concentrations
	double *minimum, //The function values
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	double *beta, //2-value array with parameters for beta prior on b2
	const int method) //ID of method to be usedfor optimization. 
{
	struct ll_param fparam = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1  b2
	unsigned type[3] =  {3,    3,  3};

	for (int i = 0; i<n_iters; i++){
		fparam.probs = probs[i];
		multimin(3, b[i], &minimum[i], type, lb, ub, &ll3f,&ll3df,&ll3fdf,
				(void *) &fparam, optim_par);
	}

}

void ll2_array_min(	
	int probs_size, //Number of vals in the prob/conc arrays
	int n_iters, //Number of different sets of probs
	double b[][2], //The vals to be optimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, //Array of concentrations
	double *minimum, //The function values
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	const int method) //ID of method to be usedfor optimization.
{
	struct ll_param fparam = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1 
	// double xmin[2]   =  {-10, 0.05};
	// double xmax[2]   =  {10,   5};
	unsigned type[2] =  {3,    3};

	for (int i = 0; i<n_iters; i++){
		fparam.probs = probs[i];
		multimin(2, b[i], &minimum[i], type, lb, ub, &ll2f,&ll2df,&ll2fdf,
				(void *) &fparam, optim_par);
	}
}

/*
This function for for calculating both the ll2 curve and the ll3 curve, and
then the AIC is used to determine whether the 2-parameter or 3-parameter curve
should be used. 
*/
void ll2_ll3_AIC(double *b, //Vector of parameters for optimization
	double *probs, //Array of probs
	double *conc, //Array of concentrations
	int probs_size, //Length of the prob/conc arrays
	double minimum, //The function value
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	double *beta, //2-value array with parameters for beta prior on b2
	const int method) //ID of method to be usedfor optimization. 
{
	struct ll_param fparam = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1  b2
	unsigned type[3] =  {3,    3,  3};

	multimin(2, b, &minimum, type, lb, ub, &ll2f, &ll2df, &ll2fdf, 
			(void *) &fparam, optim_par);
	
	double min2;
	double b2[3] = {b[0], b[1], 1.0};
	//calculate a modified log-likelihood.
	ll_all_AIC(b2, &fparam, &min2);

	multimin(3, b, &minimum, type, lb, ub, &ll3f, &ll3df, &ll3fdf, 
			(void *) &fparam, optim_par);
	//calculate a modified log-likelihood.
	double min3; 
	double b3[3] = {b[0], b[1], b[2]};
	ll_all_AIC(b3, &fparam, &min3);

	double AIC2 = 4 - min2 * 2;
	double AIC3 = 6 - min3 * 2;
	if (AIC2 < AIC3) {
		//Then 2-param is a better fit.
		for (int j = 0; j<3; j++){
			b[j] = b2[j];
		}
		minimum = min2;
	} else {
		minimum = min3;
	}

}


void array_ll2_ll3_AIC(
	int probs_size, //Number of vals in the prob/conc arrays
	int n_iters, //Number of different sets of probs
	double b[][3], //The vals to be optimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, //Array of concentrations
	double *minimum, //The function values
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	double *beta, //2-value array with parameters for beta prior on b2
	const int method) //ID of method to be usedfor optimization. 
{
	struct ll_param fparam = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1  b2
	unsigned type[3] =  {3,    3,  3};

	for (int i = 0; i<n_iters; i++){
		fparam.probs = probs[i];
		//ll2 curve fit
		multimin(2, b[i], &minimum[i], type, lb, ub, &ll2f,&ll2df,&ll2fdf,
					(void *) &fparam, optim_par);
		double min2;
		double b2[3] = {b[i][0], b[i][1], 1.0};
		ll_all_AIC(b2, &fparam, &min2);//calculate a modified log-likelihood.

		//ll3 curve fit
		multimin(3, b[i], &minimum[i], type, lb, ub, &ll3f,&ll3df,&ll3fdf,
					(void *) &fparam, optim_par);
		double min3; 
		double b3[3] = {b[i][0], b[i][1], b[i][2]};
		ll_all_AIC(b3, &fparam, &min3);//calculate a modified log-likelihood.

		//
		double AIC2 = 4 - min2 * 2;
		double AIC3 = 6 - min3 * 2;

		if (AIC2 < AIC3) {
			//Then 2-param is a better fit.
			for (int j = 0; j<3; j++){
				b[i][j] = b2[j];
			}
			minimum[i] = min2;
		} else {
			minimum[i] = min3;
		}
	}
}


void ls_driver(
	const int nparam, //number of parameters to use in dose-response curve
	const int n, //number of data points
	double *user_b, //initial guess for b (and output). Must be len 3
	double *x, //concentration values
	double *y, //probability values
	const int method) //least squares method to use
{
	gsl_vector *f = gsl_vector_alloc(n);
	gsl_vector *b = gsl_vector_alloc(nparam);
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params =
		gsl_multifit_nlinear_default_parameters();
	struct data fit_data = {.x = x,
							.y = y,
							.n = n};

	/* define LL2/LL3 least squares functions and params */
	gsl_vector_set(b, 0, user_b[0]);
	gsl_vector_set(b, 1, user_b[1]);
	fdf.n = n;
	fdf.p = nparam;
	fdf.params = &fit_data;
	if (nparam == 2){
		fdf.f = ls2_func_f;
		fdf.df = ls2_func_df;
		fdf.fvv = ls2_func_fvv;
	} else {
		gsl_vector_set(b, 2, user_b[2]);
		fdf.f = ls3_func_f;
		fdf.df = ls3_func_df;
		fdf.fvv = ls3_func_fvv;
	}

	// printf("Method: %d\n",method);
	switch(method){
		case 0:
			fdf_params.trs = gsl_multifit_nlinear_trs_lm;
			break;
		case 1:
			fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
			break;
		case 2:
			fdf_params.trs = gsl_multifit_nlinear_trs_dogleg;
			break;
		case 3:
			fdf_params.trs = gsl_multifit_nlinear_trs_ddogleg;
			break;
		case 4:
			fdf_params.trs = gsl_multifit_nlinear_trs_subspace2D;
			break;
		default:
			fdf_params.trs = gsl_multifit_nlinear_trs_lm;
	}

	solve_system(b, &fdf, &fdf_params);

	user_b[0] = gsl_vector_get(b, 0);
	user_b[1] = gsl_vector_get(b, 1);
	user_b[2] = nparam >= 3 ? gsl_vector_get(b, 2): 1.0;

	gsl_vector_free(f);
	gsl_vector_free(b);

}

void ls_array_min(
	const int nparam, //number of parameters to use in dose-response curve
	const int probs_size, //Number of vals in the prob/conc arrays
	const int n_iters, //Number of different sets of probs
	double b[][3], //The vals to be optimized. Should be of size n_iters x 3
	double *conc, //concentration values
	double probs[][probs_size], //Should be of size n_iters x probs_size
	const int method) //least squares method to use
{
	for (int i = 0; i<n_iters; i++){
		ls_driver(nparam, probs_size, b[i], conc, probs[i], method);
	}

}

void cls3_min(double *b, //Vector of parameters for optimization
	double *probs, //Array of probs
	double *conc, //Array of concentrations
	int probs_size, //Length of the prob/conc arrays
	double minimum, //The function value
	double *lb, //Lower bound for constrained optim
	double *ub) //Upper bound for constrained optim
{

	struct ll_param fparam = { .probs = probs, 
								.conc = conc, 
								.probs_size = probs_size
							};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,6,0};

	//                   b0	  b1  b2
	unsigned type[3] =  {3,    3,  3};

	multimin(3, b, &minimum, type, lb, ub, &cls3f_error,NULL,NULL,
			(void *) &fparam, optim_par);
}


void cls3_array_min(
	int probs_size, //Number of vals in the prob/conc arrays
	int n_iters, //Number of different sets of probs
	double b[][3], //The vals to be optimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, //Array of concentrations
	double *minimum, //The function values
	double *lb, //Lower bound for constrained optim
	double *ub) //Upper bound for constrained optim
{
	struct ll_param fparam ={ .probs = probs[0], 
							  .conc = conc, 
							  .probs_size = probs_size
							};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,6,0};

	//                   b0	  b1  b2
	unsigned type[3] =  {3,    3,  3};

	for (int i = 0; i<n_iters; i++){
		fparam.probs = probs[i];
		multimin(3, b[i], &minimum[i], type, lb, ub, &cls3f_error,NULL,NULL,
				(void *) &fparam, optim_par);
	}

}