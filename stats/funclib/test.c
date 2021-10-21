#include "multimin.h"
#include "ll_mod.h"

int main(){

	printf("Starting!\n");

	double conc[21] = {8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125};
	double probs[21] = {0, 0, 0, 0, 0.6875, 0.846153846153846, 1, 0.857142857142857, 0.923076923076923, 0.954545454545455, 0.875, 0, 0, 0, 0, 0, 0, 0, 0.083333333333333, 0.15, 0.916666666666667};
	double b[3] = {.90,-1.191,.92};

	struct ll3_param fparam = { .probs = probs, .conc = conc, .sigsquare = 1e6, .probs_size = 21, .weib = {2.0,1.0}};

	double fval;	

	double grad[3] = {0,0,0};

	ll3f(3, b, &fparam, &fval);

	ll3df(3, b, &fparam, grad);

	printf("%f\n", fval);

	printf("%f, %f, %f\n", grad[0], grad[1], grad[2]);

	return 0;
}