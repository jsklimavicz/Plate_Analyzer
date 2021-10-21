/* ll.c
 * Contains the log-likelihood functions and jacobians for curve-fitting.
 * Also contains the error function for least squares fitting. 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "ll_mod.h"

/*
GLOBALS
*/

double LB = 1.E-10;
double DEF_SS = 1.0E6;
double DEF_BETA[2] = {1.5,1.01};
double DEF_RETURN = 1E10;

#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

void ll3f(const size_t n, 
	const double *b, 
	void * fparams, 
	double * fval) {
	/*
	Log-likelihood function of the three parameter dose-response curve 
						    b2
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
	b2 ~ Beta(beta_params).
	*/

	struct ll3_param *par = (struct ll3_param *) fparams;

	const double * probs = par->probs;
	const double * conc = par->conc;
	const double sigsquare = par->sigsquare;
	const int probs_size = par->probs_size;
	const double * beta = par->beta;


	if (b[2] < LB) { 
		*fval = DEF_RETURN;
		return; 
	}

	double alpha;

	double lb2 = log(b[2]);

	*fval = (b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare) //MVN prior
		- (beta[0]-1)*lb2 - (beta[1]-1)*log(1-b[2]); //Beta prior
	// printf("Before: ll = %f \n", *fval);
	// Only iterate through probs and conc arrays once and
	// add to log-likelihood each time. 
	for (int i = 0; i < probs_size; i++) {
		alpha = 1.0 + pow(M_E, b[0] + b[1]*conc[i]);
		if (alpha - b[2] < LB) { 
			*fval = DEF_RETURN;
			return; 
		}
		*fval -= probs[i]*(log(alpha-b[2]) - lb2) - log(alpha);
		// printf("After %i: ll = %f \n", i, *fval);
	}
	*fval -= lb2 * probs_size;
	// printf(" %f \n", *fval);
}


void ll3df(const size_t n, 
	const double * b,
	void * fparams,
	double * grad) 
{
	/*
	Jacobian of the log-likelihood function  of the three parameter 
	dose-response curve 
						    b2
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein prior is b0, b1 ~ MVN(0, sigma*I2) and 
	b2 ~ Beta(beta_params).
	*/
	struct ll3_param *par = (struct ll3_param *) fparams;

	const double * probs = par->probs;
	const double * conc = par->conc;
	const double sigsquare = par->sigsquare;
	const int probs_size = par->probs_size;
	const double * beta = par->beta;

	grad[0] = b[0]/sigsquare;
	grad[1] = b[1]/sigsquare;
	grad[2] = -(beta[0]-1)/b[2] + (beta[1]-1)/(1-b[2]);

	double xi, alpha, d, m;
	// Only iterate through probs and conc arrays once and
	// add to each component of the gradient each time. 
	for (int i = 0; i < probs_size; i++) {
		xi = pow(M_E, b[0] + b[1]*conc[i]);
		alpha = 1.0 + xi;
		d = alpha - b[2];
		m = probs[i]*xi / d - xi / alpha;
		grad[0] -= m;
		grad[1] -= conc[i]*m;
		grad[2] -= (1-probs[i])/b[2] - probs[i]/d;
	}
}

void ll3fdf(const size_t n, 
	const double *b,
	void * fparams,
	double *fval,
	double *grad) 
{
	// ll3f(n, b, fparams, fval);
	// ll3df(n, b, fparams, grad);

	struct ll3_param *par = (struct ll3_param *) fparams;

	const double * probs = par->probs;
	const double * conc = par->conc;
	const double sigsquare = par->sigsquare;
	const int probs_size = par->probs_size;
	const double * beta = par->beta;

	if (b[2] < LB) { 
		*fval = DEF_RETURN;
		return; 
	}

	double lb2 = log(b[2]);

	*fval = (b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare) //MVN prior
		- (beta[0]-1)*lb2 - (beta[1]-1)*log(1-b[2]); //Beta prior

	grad[0] = b[0]/sigsquare;
	grad[1] = b[1]/sigsquare;
	grad[2] = -(beta[0]-1)/b[2] + (beta[1]-1)/(1-b[2]);

	double xi, alpha, d, m;
	// Only iterate through probs and conc arrays once and
	// add to each component of the gradient each time. 
	for (int i = 0; i < probs_size; i++) {
		xi = pow(M_E, b[0] + b[1]*conc[i]);
		alpha = 1.0 + xi;
		if (alpha - b[2] < LB) { 
			*fval = DEF_RETURN;
		}
		d = alpha - b[2];
		m = probs[i]*xi / d - xi / alpha;
		grad[0] -= m;
		grad[1] -= conc[i]*m;
		grad[2] -= (1-probs[i])/b[2] - probs[i]/d;
		*fval -= probs[i]*(log(d) - lb2) - log(alpha);
	}
	*fval -= lb2 * probs_size;
}



void ll2f(const size_t n, 
	const double *b, 
	void * fparams, 
	double * fval) {
	/*
	Log-likelihood function of the two parameter dose-response curve 
						    1
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein priors are b0, b1 ~ MVN(0, sigma*I2) 
	*/

	struct ll2_param *par = (struct ll2_param *) fparams;

	const double * probs = par->probs;
	const double * conc = par->conc;
	const double sigsquare = par->sigsquare;
	const int probs_size = par->probs_size;

	*fval = (b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare); //MVN prior

	double xi, x;
	// Only iterate through probs and conc arrays once and
	// add to each component of the gradient each time. 
	for (int i = 0; i < probs_size; i++) {
		x = b[0] + b[1]*conc[i];
		xi = pow(M_E, x);
		*fval -= x*probs[i] - log(1.0 + xi);
	}

}


void ll2df(const size_t n, 
	const double * b,
	void * fparams,
	double * grad) 
{
	/*
	Jacobian of the log-likelihood function  of the two parameter 
	dose-response curve 
						    b2
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein prior is b0, b1 ~ MVN(0, sigma*I2).
	*/
	struct ll2_param *par = (struct ll2_param *) fparams;

	const double * probs = par->probs;
	const double * conc = par->conc;
	const double sigsquare = par->sigsquare;
	const int probs_size = par->probs_size;

	grad[0] = b[0]/sigsquare;
	grad[1] = b[1]/sigsquare;

	double xi, d;
	// Only iterate through probs and conc arrays once and
	// add to each component of the gradient each time. 
	for (int i = 0; i < probs_size; i++) {
		xi = pow(M_E, b[0] + b[1]*conc[i]);
		d = probs[i] - xi/(1.0 + xi);
		grad[0] -= d ;
		grad[1] -= conc[i]*d;
	}
}

void ll2fdf(const size_t n, 
	const double *b,
	void * fparams,
	double *fval,
	double *grad) 
{
	// ll3f(n, b, fparams, fval);
	// ll3df(n, b, fparams, grad);

	struct ll2_param *par = (struct ll2_param *) fparams;

	const double * probs = par->probs;
	const double * conc = par->conc;
	const double sigsquare = par->sigsquare;
	const int probs_size = par->probs_size;

	*fval = (b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare); //MVN prior

	grad[0] = b[0]/sigsquare;
	grad[1] = b[1]/sigsquare;

	double xi, alpha, x, d;
	// Only iterate through probs and conc arrays once and
	// add to each component of the gradient each time. 
	for (int i = 0; i < probs_size; i++) {
		x = b[0] + b[1]*conc[i];
		xi = pow(M_E, x);
		alpha = 1.0 + xi;
		d = probs[i] - xi/alpha;

		grad[0] -= d ;
		grad[1] -= conc[i]*d;

		*fval -= x*probs[i] - log(alpha);
	}

}
