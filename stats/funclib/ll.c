#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* globals */
double LB = 1.E-10;
double DEF_SS = 1.0E6;
double DEF_WEIB[2] = {2.0, 1.0};
double DEF_RETURN = 1E10;

double ll3_full(const double *b, 
	const double *probs, const double *conc,
	const double sigsquare, const double *weib)
{
	if (b[2] < LB) { return DEF_RETURN; }

	double alpha;
	double wk = weib[0];
	double wl = weib[1] * pow(wk/(wk-1.0), 1/wk);
	double lb2 = log(b[2]);

	double ll = -(b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare)
		- pow((b[2]/wl),wk) + (wk-1)*log(b[2]);

	int size = sizeof(probs)/sizeof(probs[0]);
	for (int i = 0; i < size; i++) {
		alpha = 1.0 + pow(M_E, b[0] + b[1]*conc[i]);
		if (alpha - b[2] < LB) { return 1E10; }
		ll += probs[i]*(log(alpha-b[2]) - lb2) - log(alpha);
	}

	ll += lb2 * size;
	return -ll;
}


double ll2_full(const double *b, 
	const double *probs, const double *conc,
	const double sigsquare)
{

	double alpha;
	double ll = -(b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare);
	int size = sizeof(probs)/sizeof(probs[0]);
	for (int i = 0; i < size; i++) {
		alpha = log(1.0 + pow(M_E, b[0] + b[1]*conc[i]));
		ll += (b[0] + b[1] * conc[i]) * probs[i];
	}

	return -ll;
}


double *ll2_jac_full(const double *b, 
	 const double *probs, const double *conc,
	 const double sigsquare) 
{
	double *g = malloc(2*sizeof(double));

	g[0] = b[0]/sigsquare;
	g[1] = b[1]/sigsquare;

	int size = sizeof(probs)/sizeof(probs[0]);

	double alpha;
	for (int i = 0; i < size; i++) {
		alpha = pow(M_E, b[0] + b[1]*conc[i]);
		alpha = alpha/(1.0 + alpha);
		g[0] -= probs[i] - alpha;
		g[1] -= conc[i]*probs[i] - alpha * conc[i];
	}
	return g;
}

double ls(const double *b,const double *probs, 
	const double *conc) 
{
	bsize = sizeof(b)/sizeof(b[0]);
}






double ll3(const double *b, 
	 const double *probs,
	 const double *conc) 
{
	return ll3_full(b, probs, conc, DEF_SS, DEF_WEIB);
}

double ll2(const double *b, 
	 const double *probs,
	 const double *conc) 
{
	return ll2_full(b, probs, conc, DEF_SS);
}
