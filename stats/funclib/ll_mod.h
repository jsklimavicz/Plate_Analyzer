/*
 * Contains the log-likelihood functions and jacobians for curve-fitting.
 * Also contains the error function for least squares fitting. 
 */



struct ll2_param {
  double *probs;
  double *conc;
  int probs_size;
  double sigsquare;
};

struct ll3_param {
  double *probs;
  double *conc;
  int probs_size;
  double sigsquare;
  double *beta ;
};

struct ls_param {
  double *probs;
  double *conc;
  int probs_size;
};


/*Three-parameter log-likelihood*/
void ll3f(const size_t n, 
	const double *x, 
	void * fparams, 
	double *fval);

void ll3df(const size_t n, 
  const double * b,
  void * fparams,
  double * grad) ;

void ll3fdf(const size_t n, 
  const double * b,
  void * fparams,
  double *fval,
  double *grad) ;


/*Two-parameter log-likelihood*/
void ll2f(const size_t n, 
	const double *x, 
	void *fparams, 
	double *fval);

void ll2df(const size_t n, 
	const double *x,
	void *fparams,
	double *grad) ;

void ll2fdf(const size_t n, 
	const double *x,
	void *fparams,
	double *fval,
	double *grad) ;


// double ls(const double *b,
// 	const double *probs, 
// 	const double *conc, 
// 	const int probs_size, 
// 	const int num_param) ;

// double ls2(const double *b,const double *probs, 
// 	const double *conc, const int probs_size) ;

// double ls3(const double *b,const double *probs, 
// 	const double *conc, const int probs_size) ;
