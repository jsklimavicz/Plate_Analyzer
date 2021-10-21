/* ll.c
 * Contains the log-likelihood functions and jacobians for curve-fitting.
 * Also contains the error function for least squares fitting. 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "ll.h"

double ll3_full(const double *b, 
	const double *probs, 
	const double *conc, 
	const int probs_size, 
	const double sigsquare, 
	const double *weib)
{
	/*
	Log-likelihood function of the three parameter dose-response curve 
						    b2
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
	b2 ~ Weibull(weibull_params).
	*/
	if (b[2] < LB) { return DEF_RETURN; }

	double alpha;

	double wk = weib[0];
	double wl = weib[1] * pow(wk/(wk-1.0), 1/wk);
	double lb2 = log(b[2]);
	double ll = (b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare)
		+ pow((b[2]/wl),wk) + (wk-1)*log(b[2]);

	// printf("Before start: ll = %f \n", ll);

	// Only iterate through probs and conc arrays once and
	// add to log-likelihood each time. 
	for (int i = 0; i < probs_size; i++) {
		alpha = 1.0 + pow(M_E, b[0] + b[1]*conc[i]);
		if (alpha - b[2] < LB) { return 1E10; }
		ll -= probs[i]*(log(alpha-b[2]) - lb2) - log(alpha);
		ll -= lb2;
		// printf("After %i: ll = %f \n", i, ll);
	}

	// ll += lb2 * size;

	return ll;
}


double ll2_full(const double *b, 
	const double *probs, 
	const double *conc, 
	const int probs_size, 
	const double sigsquare)
{
	/*
	Log-likelihood function of the three parameter dose-response curve 
						    1
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein priors are b0, b1 ~ MVN(0, sigma*I2).
	*/
	double alpha;
	double ll = (b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare);
	// Only iterate through probs and conc arrays once and
	// add to log-likelihood each time. 
	for (int i = 0; i < probs_size; i++) {
		alpha = log(1.0 + pow(M_E, b[0] + b[1]*conc[i]));
		ll -= (b[0] + b[1] * conc[i]) * probs[i] - alpha;
	}
	return ll;
}


void ll2_jac_full(const double *b, 
	 const double *probs, 
	 const double *conc, 
	 const int probs_size, 
	 double *g, //returned vector
	 const double sigsquare) 
{
	/*
	Jacobian of the log-likelihood function  of the two parameter 
	dose-response curve 
						    1
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein prior is b0, b1 ~ MVN(0, sigma*I2).
	*/

	g[0] = b[0]/sigsquare;
	g[1] = b[1]/sigsquare;

	double alpha;
	// Only iterate through probs and conc arrays once and
	// add to each component of the gradient each time. 
	for (int i = 0; i < probs_size; i++) {
		alpha = pow(M_E, b[0] + b[1]*conc[i]);
		alpha = alpha/(1.0 + alpha);
		g[0] -= probs[i] - alpha;
		g[1] -= conc[i]*probs[i] - alpha * conc[i];
	}
}

void ll3_jac_full(const double *b, 
	 const double *probs, 
	 const double *conc, 
	 const int probs_size, 
	 double *g, //returned vector
	 const double sigsquare, 
	 const double *weib) 
{
	/*
	Jacobian of the log-likelihood function  of the two parameter 
	dose-response curve 
						    b2
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein prior is b0, b1 ~ MVN(0, sigma*I2) and 
	b2 ~ Weibull(weibull_params).
	*/
	double wk = weib[0];
	double wl = weib[1] * pow(wk/(wk-1.0), 1/wk);

	// double *g = malloc(3*sizeof(double));

	g[0] = b[0]/sigsquare;
	g[1] = b[1]/sigsquare;
	g[2] = -(wk-1)/b[2] + (wk/b[2]) * pow((b[2]/wl),wk);

	double xi;
	double alpha;
	double d;
	double m; 
	double l;
	// Only iterate through probs and conc arrays once and
	// add to each component of the gradient each time. 
	for (int i = 0; i < probs_size; i++) {
		xi = pow(M_E, b[0] + b[1]*conc[i]);
		alpha = 1.0 + xi;
		d = alpha - b[2];
		m = probs[i]*xi / d;
		l = xi / alpha;
		g[0] -= m - l;
		g[1] -= conc[i]*(m - l);
		g[2] -= (1-probs[i])/b[2] - probs[i]/d;
	}
}

double ls(const double *b,
	const double *probs, 
	const double *conc, 
	const int probs_size, 
	const int num_param) 
{
	/*
	Dose-response curve for least-squares fitting. If len(b)=2, then 
	b2 = 1; else if len(b)=3, then b2 = b[2]. Returns error
					  |	    	b2				| 2
			  error = | ------------------ - y 	|
					  | 1 + exp(b0 + b1*x)		|
	*/
	double error = 0.0;
	double temp = 0.0;
	double b2 = (probs_size == 2 ? 1.0 : b[2] );

	for (int i = 0; i < probs_size; i++) {
		temp = b2/(1.0 + pow(M_E, b[0] + b[1]*conc[i])) - probs[i];
		error += temp * temp;
	}
	return error;
}


/* 
WRAPPER FUNCTIONS
*/
double ll3(const double *b, 
	 const double *probs,
	 const double *conc, 
	 const int probs_size) 
{
	return ll3_full(b, probs, conc, probs_size, DEF_SS, DEF_WEIB);
}

double ll2(const double *b, 
	 const double *probs,
	 const double *conc, 
	 const int probs_size) 
{
	return ll2_full(b, probs, conc, probs_size, DEF_SS);
}

void ll2j(const double *b, 
	 const double *probs, 
	 const double *conc, 
	 const int probs_size, 
	 double *g) //g is the returned gradient.
{
	ll2_jac_full(b, probs, conc, probs_size, g, DEF_SS);
}

void ll3j(const double *b, 
	 const double *probs, 
	 const double *conc, 
	 const int probs_size, 
	 double *g) //g is the returned gradient.
{

	ll3_jac_full(b, probs, conc, probs_size, g, DEF_SS, DEF_WEIB);
}

double ls2(const double *b,const double *probs, 
	const double *conc, const int probs_size) 
{
	return ls(b, probs, conc, probs_size, 2);
}

double ls3(const double *b,const double *probs, 
	const double *conc, const int probs_size) 
{
	return ls(b, probs, conc, probs_size, 3);
}
