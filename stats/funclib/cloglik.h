


void ll3_min(double *b, //The minimal value that is found.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare, 
	double *beta);

void ll2_min(double *b, //The minimal value that is found.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare);

void ll3_array_min(
	int probs_size, 
	int n_iters,
	double b[][n_iters], //The minima value that are found. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double minimum, //The function value
	double sigsquare, 
	double *beta);

void ll2_array_min(
	int probs_size, 
	int n_iters,
	double b[][n_iters], //The minima value that are found. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double minimum, //The function value
	double sigsquare);