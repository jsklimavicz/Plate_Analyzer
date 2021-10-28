/*
cloglik.c
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

void ll3_min(double *b, //The minimal value that is found.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare, 
	double *beta)
{
	struct ll3_param fparam = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
  	struct multimin_params optim_par = {.1,1e-2,250,1e-3,1e-5,5,0};

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
	double sigsquare)
{
	struct ll2_param fparam = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size};
  	struct multimin_params optim_par = {.1,1e-2,20,1e-3,1e-5,5,0};
	
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
	double *beta)
{
	// printf("%i\n", probs_size);
	// printf("%i\n", n_iters);
	struct ll3_param fparam = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
  	struct multimin_params optim_par = {.1,1e-2,250,1e-3,1e-5,5,0};

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
	double sigsquare)
{
	struct ll2_param fparam = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size};
  	struct multimin_params optim_par = {.1,1e-2,20,1e-3,1e-5,5,0};

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
	double *beta)
{
	struct ll3_param fparam3 = { .probs = probs, 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
	struct multimin_params optim_par = {.1,1e-2,250,1e-3,1e-5,5,0};

	//                   b0	  b1  b2
	double xmin[3]   =  {-10, -5,  0};
	double xmax[3]   =  {10,   5,  1};
	unsigned type[3] =  {3,    3,  3};

	multimin(2, b, &minimum, type, xmin, xmax, &ll2f, &ll2df, &ll2fdf, 
			(void *) &fparam3, optim_par);
	
	double min2;
	double b2[3] = {b[0], b[1], 1.0};
	//calculate a modified log-likelihood.
	ll_all_AIC(b2, &fparam3, &min2);

	multimin(3, b, &minimum, type, xmin, xmax, &ll3f, &ll3df, &ll3fdf, 
			(void *) &fparam3, optim_par);
	//calculate a modified log-likelihood.
	double min3; 
	double b3[3] = {b[0], b[1], b[2]};
	ll_all_AIC(b3, &fparam3, &min3);

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
	double *beta)
{
	struct ll3_param fparam3 = { .probs = probs[0], 
								.conc = conc, 
								.sigsquare = sigsquare, 
								.probs_size = probs_size, 
								.beta = beta};
	struct multimin_params optim_par = {.1,1e-2,250,1e-3,1e-5,5,0};

	//                   b0	  b1  b2
	double xmin[3]   =  {-10, -5,  0};
	double xmax[3]   =  {10,   5,  1};
	unsigned type[3] =  {3,    3,  3};

	for (int i = 0; i<n_iters; i++){
		fparam3.probs = probs[i];
		//ll2 curve fit
		multimin(2, b[i], &minimum[i], type, xmin, xmax,&ll2f,&ll2df,&ll2fdf,
					(void *) &fparam3, optim_par);
		double min2;
		double b2[3] = {b[i][0], b[i][1], 1.0};
		ll_all_AIC(b2, &fparam3, &min2);//calculate a modified log-likelihood.

		//ll3 curve fit
		multimin(3, b[i], &minimum[i], type, xmin, xmax, &ll3f,&ll3df,&ll3fdf,
					(void *) &fparam3, optim_par);
		double min3; 
		double b3[3] = {b[i][0], b[i][1], b[i][2]};
		ll_all_AIC(b3, &fparam3, &min3);//calculate a modified log-likelihood.

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