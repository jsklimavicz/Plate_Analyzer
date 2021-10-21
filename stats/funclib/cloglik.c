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
	struct ll3_param fparam = { .probs = probs, .conc = conc, .sigsquare = sigsquare, .probs_size = probs_size, .beta = beta};
  	struct multimin_params optim_par = {.1,1e-2,250,1e-3,1e-5,5,0};

	double xmin[3], xmax[3];
	unsigned type[3];

	type[0]=3;
	xmin[0]=10;
	xmax[0]=-10;

	type[1]=3;
	xmin[1]=-5;
	xmax[1]=5;

	type[2]=3;
	xmin[2]=0;
	xmax[2]=1;

	multimin(3, b, &minimum, type, xmin, xmax, &ll3f,&ll3df,&ll3fdf, (void *) &fparam, optim_par);
}

void ll2_min(double *b, //The minimal value that is found.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare)
{
	struct ll2_param fparam = { .probs = probs, .conc = conc, .sigsquare = sigsquare, .probs_size = probs_size};
  	struct multimin_params optim_par = {.1,1e-2,20,1e-3,1e-5,5,10};

	double xmin[2], xmax[2];
	unsigned type[2];

	type[0]=3;
	xmin[0]=10;
	xmax[0]=-10;

	type[1]=3;
	xmin[1]=-5;
	xmax[1]=5;

	multimin(2, b, &minimum, type, xmin, xmax,&ll3f,&ll3df,&ll3fdf, (void *) &fparam, optim_par);
}

void ll3_array_min(
	int probs_size, 
	int n_iters,
	double b[][n_iters], //The minima value that are found. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double minimum, //The function value
	double sigsquare, 
	double *beta)
{
	// printf("%i\n", probs_size);
	// printf("%i\n", n_iters);
	struct ll3_param fparam = { .probs = probs[0], .conc = conc, .sigsquare = sigsquare, .probs_size = probs_size, .beta = beta};
  	struct multimin_params optim_par = {.1,1e-2,250,1e-3,1e-5,5,0};

	double xmin[3], xmax[3];
	unsigned type[3];

	type[0]=3;
	xmin[0]=10;
	xmax[0]=-10;

	type[1]=3;
	xmin[1]=-5;
	xmax[1]=5;

	type[2]=3;
	xmin[2]=0;
	xmax[2]=1;

	for (int i = 0; i<n_iters; i++){
		// for(int j =0; j<probs_size; j++){
		// 	printf("%f", probs[i][j]);
		// }
		// printf("\n");
		fparam.probs = probs[i];
		multimin(3, b[i], &minimum, type, xmin, xmax, &ll3f,&ll3df,&ll3fdf, (void *) &fparam, optim_par);
		// for(int j =0; j<3; j++){
		// 	printf("%f ", b[i][j]);
		// }
		// printf("\n");
	}

}

void ll2_array_min(	
	int probs_size, 
	int n_iters,
	double b[][n_iters], //The minima value that are found. Should be of size n_iters x 3
	double probs[][n_iters], //Should be of size n_iters x probs_size
	double *conc, 
	double minimum, //The function value
	double sigsquare)
{
	struct ll2_param fparam = { .probs = probs[0], .conc = conc, .sigsquare = sigsquare, .probs_size = probs_size};
  	struct multimin_params optim_par = {.1,1e-2,20,1e-3,1e-5,5,10};

	double xmin[2], xmax[2];
	unsigned type[2];

	type[0]=3;
	xmin[0]=10;
	xmax[0]=-10;

	type[1]=3;
	xmin[1]=-5;
	xmax[1]=5;

	for (int i = 0; i<n_iters; i++){
		fparam.probs = probs[i];
		multimin(2, b[i], &minimum, type, xmin, xmax, &ll3f,&ll3df,&ll3fdf, (void *) &fparam, optim_par);
	}
}