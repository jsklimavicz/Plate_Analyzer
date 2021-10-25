#include "multimin.h"
#include "ll_mod.h"

int main(){

	double conc[20] = {8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125};
	double probs[20] = {0, 0, 0, .33, 0.6875, 0.846153846153846, 1, 0.857142857142857, 0.923076923076923, 0.914545454545455, 0, 0, 0, 0, 0, .2, .43, 0.583333333333333, 0.85, 0.916666666666667};
	double b[3] = {.89,-1.3,.917};

	struct ll3_param fparam = { .probs = probs, .conc = conc, .sigsquare = 1e6, .probs_size = 20, .beta = {1.5,1.01}};

	double fval;	

	double grad[3] = {0,0,0};

	ll3f(3, b, &fparam, &fval);
	printf("%f\n", fval);

	ll3df(3, b, &fparam, grad);

	printf("%f, %f, %f\n", grad[0], grad[1], grad[2]);

	ll3fdf(3, b, &fparam, &fval, grad);
	printf("%f\n", fval);
	printf("%f, %f, %f\n", grad[0], grad[1], grad[2]);

  	double minimum;
  	struct multimin_params optim_par = {.1,1e-2,20,1e-3,1e-5,5,10};

	  double xmin[3], xmax[3];
	  unsigned type[3];

  	  /* minimum constrained in [-1,2]x(-5,5] */
	type[0]=3;
	xmin[0]=10;
	xmax[0]=-10;

	type[1]=3;
	xmin[1]=-5;
	xmax[1]=5;

	type[2]=3;
	xmin[2]=0;
	xmax[2]=1;

  	optim_par.verbosity=10;

  	printf("%u", optim_par.method);

  	// multimin(3, b, &minimum,type,xmin,xmax,&ll3f,&ll3df,&ll3fdf, (void *) &fparam, optim_par);
  	multimin(3, b, &minimum,NULL,NULL,NULL,&ll3f,&ll3df,&ll3fdf, (void *) &fparam, optim_par);

  	printf("unconstrained minimum %e at [%f, %f, %f]\n",minimum,b[0],b[1],b[2]);

	return 0;
}