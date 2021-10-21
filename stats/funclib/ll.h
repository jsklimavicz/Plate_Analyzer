/*
 * Contains the log-likelihood functions and jacobians for curve-fitting.
 * Also contains the error function for least squares fitting. 
 */

/*
GLOBALS
*/
double LB = 1.E-10;
double DEF_SS = 1.0E6;
double DEF_WEIB[2] = {2.0, 1.0};
double DEF_RETURN = 1E10;

#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

double ll3_full(const double *b, 
	const double *probs, 
	const double *conc, 
	const int probs_size, 
	const double sigsquare, 
	const double *weib);

double ll2_full(const double *b, 
	const double *probs, 
	const double *conc, 
	const int probs_size, 
	const double sigsquare);

void ll2_jac_full(const double *b, 
	 const double *probs, 
	 const double *conc, 
	 const int probs_size, 
	 double *g, //returned vector
	 const double sigsquare) ;

void ll3_jac_full(const double *b, 
	 const double *probs, 
	 const double *conc, 
	 const int probs_size, 
	 double *g, //returned vector
	 const double sigsquare, 
	 const double *weib) ;

double ls(const double *b,
	const double *probs, 
	const double *conc, 
	const int probs_size, 
	const int num_param) ;

double ll3(const double *b, 
	 const double *probs,
	 const double *conc, 
	 const int probs_size) ;

double ll2(const double *b, 
	 const double *probs,
	 const double *conc, 
	 const int probs_size) ;

void ll2j(const double *b, 
	 const double *probs, 
	 const double *conc, 
	 const int probs_size, 
	 double *g);

double ls2(const double *b,const double *probs, 
	const double *conc, const int probs_size) ;

double ls3(const double *b,const double *probs, 
	const double *conc, const int probs_size) ;
