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

void ll3_min(double *b, //The minimal value that is found.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare, 
	double *beta, 
	const int method)
{
	struct ll_param fparam = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1  b2
	double xmin[3]   =  {-10, -5,  0};
	double xmax[3]   =  {10,   5,  1};
	unsigned type[3] =  {3,    3,  3};

	multimin(3, b, &minimum, type, xmin, xmax, &ll3f,&ll3df,&ll3fdf,
			(void *) &fparam, optim_par);
}

void ll2_min(double *b, //The minimal value that is found.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare,
	const int method)
{
	struct ll_param fparam = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};
	
	//                   b0	  b1 
	double xmin[2]   =  {-10, -5};
	double xmax[2]   =  {10,   5};
	unsigned type[2] =  {3,    3};

	multimin(2, b, &minimum, type, xmin, xmax,&ll2f,&ll2df,&ll2fdf,
			(void *) &fparam, optim_par);
}

void ll3_array_min(
	int probs_size, 
	int n_iters,
	double b[][3], //The val to be minimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double *minimum, //The function value
	double sigsquare, 
	double *beta,
	const int method)
{
	// printf("%i\n", probs_size);
	// printf("%i\n", n_iters);
	struct ll_param fparam = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1  b2
	double xmin[3]   =  {-10, -5,  0};
	double xmax[3]   =  {10,   5,  1};
	unsigned type[3] =  {3,    3,  3};

	for (int i = 0; i<n_iters; i++){
		fparam.probs = probs[i];
		// for (int j = 0; j < 3; j++){
		// 	printf("%.5f ", b[i][j]);
		// }
		// printf("\n");
		multimin(3, b[i], &minimum[i], type, xmin, xmax, &ll3f,&ll3df,&ll3fdf,
				(void *) &fparam, optim_par);
	}

}

void ll2_array_min(	
	int probs_size, 
	int n_iters,
	double b[][2], //The vals to be minimized. Should be of size n_iters x 2
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double *minimum, //The function value
	double sigsquare,
	const int method)
{
	struct ll_param fparam = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size};
  	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1 
	double xmin[2]   =  {-10, -5};
	double xmax[2]   =  {10,   5};
	unsigned type[2] =  {3,    3};

	for (int i = 0; i<n_iters; i++){
		fparam.probs = probs[i];
		multimin(2, b[i], &minimum[i], type, xmin, xmax, &ll2f,&ll2df,&ll2fdf,
				(void *) &fparam, optim_par);
	}
}


void ll2_ll3_AIC(double *b, //The val to be minimized. Must have length 3.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare,
	double *beta,
	const int method)
{
	struct ll_param fparam = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1  b2
	double xmin[3]   =  {-10, -5,  0};
	double xmax[3]   =  {10,   5,  1};
	unsigned type[3] =  {3,    3,  3};

	multimin(2, b, &minimum, type, xmin, xmax, &ll2f, &ll2df, &ll2fdf, 
			(void *) &fparam, optim_par);
	
	double min2;
	double b2[3] = {b[0], b[1], 1.0};
	//calculate a modified log-likelihood.
	ll_all_AIC(b2, &fparam, &min2);

	multimin(3, b, &minimum, type, xmin, xmax, &ll3f, &ll3df, &ll3fdf, 
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
	int probs_size, 
	int n_iters,
	double b[][3], //The vals to be minimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double *minimum, //The function value
	double sigsquare, 
	double *beta,
	const int method)
{
	struct ll_param fparam = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
	struct multimin_params optim_par = {.1,1e-2,500,1e-3,1e-5,method,0};

	//                   b0	  b1  b2
	double xmin[3]   =  {-10, -5,  0};
	double xmax[3]   =  {10,   5,  1};
	unsigned type[3] =  {3,    3,  3};

	for (int i = 0; i<n_iters; i++){
		fparam.probs = probs[i];
		//ll2 curve fit
		multimin(2, b[i], &minimum[i], type, xmin, xmax,&ll2f,&ll2df,&ll2fdf,
					(void *) &fparam, optim_par);
		double min2;
		double b2[3] = {b[i][0], b[i][1], 1.0};
		ll_all_AIC(b2, &fparam, &min2);//calculate a modified log-likelihood.

		//ll3 curve fit
		multimin(3, b[i], &minimum[i], type, xmin, xmax, &ll3f,&ll3df,&ll3fdf,
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


void ls_driver(const int nparam, //number of parameters to use
	const int n, //number of data points
	double *user_b, //initial guess for b (and output). Must be len 3
	double *x, //concentration values
	double *y, // y-values
	const int method) //for picking solving method
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
	const int nparam,
	const int probs_size, 
	const int n_iters,
	double b[][3], //The vals to be minimized. Should be of size n_iters x 3
	double *conc,
	double probs[][probs_size], //Should be of size n_iters x probs_size
	const int method)
{
	for (int i = 0; i<n_iters; i++){

		ls_driver(nparam, probs_size, b[i], conc, probs[i], method);
	}

}

