#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* globals */
double LB = 1.E-10;
double DEF_SS = 1.0E6;
double DEF_WEIB[2] = {2.0, 1.0};
double DEF_RETURN = 1E10;

double ll3_full(const int n_prbs, const double *b, 
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

	for (int i = 0; i < n_prbs; i++) {
		alpha = 1.0 + pow(M_E, b[0] + b[1]*conc[i]);
		if (alpha - b[2] < LB) { return 1E10; }
		ll += probs[i]*(log(alpha-b[2]) - lb2) - log(alpha);
	}

	ll += lb2 * n_prbs;
	return -ll;
}


double ll2_full(const int n_prbs, const double *b, 
	const double *probs, const double *conc,
	const double sigsquare)
{

	double alpha;
	double ll = -(b[0] * b[0] + b[1]*b[1])/(2.0*sigsquare);

	for (int i = 0; i < n_prbs; i++) {
		alpha = log(1.0 + pow(M_E, b[0] + b[1]*conc[i]));
		ll += (b[0] + b[1] * conc[i]) * probs[i];
	}

	return -ll;
}


double *ll2_jac_full(const int n_prbs, const double *b, 
	 const double *probs, const double *conc,
	 const double sigsquare) 
{
	double *g = malloc(2*sizeof(double));

	g[0] = b[0]/sigsquare;
	g[1] = b[1]/sigsquare;

	double alpha;
	for (int i = 0; i < n_prbs; i++) {
		alpha = pow(M_E, b[0] + b[1]*conc[i]);
		alpha = alpha/(1.0 + alpha);
		g[0] -= probs[i] - alpha;
		g[1] -= conc[i]*probs[i] - alpha * conc[i];
	}
	return g;
}








double ll3(const int n_prbs,
	 const double *b, 
	 const double *probs,
	 const double *conc) 
{
	return ll3_full(n_prbs, b, probs, conc, DEF_SS, DEF_WEIB);
}

double ll2(const int n_prbs,
	 const double *b, 
	 const double *probs,
	 const double *conc) 
{
	return ll2_full(n_prbs, b, probs, conc, DEF_SS);
}

// double* ll2_jac(const int n_prbs,
// 	 const double *b, 
// 	 const double *probs,
// 	 const double *conc) 
// {
// 	return ll2_full(n_prbs, b, probs, conc, DEF_SS);
// }