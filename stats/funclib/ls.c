/*
stats/funclib/ls.c
This file is a component of Plate_Analyzer, which can count Drosophila
L1 larvae, classify them as alive or dead, and determine dose-repsonse
information based on live/dead count data. 

Copyright (C) 2021 James Klimavicz

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include "ls.h"
#include "ll_mod.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>


/* 2-param Log-logit curve: 1 / (1 + exp(b0 + b1 * x) ) */
double ls2_fun(const double b0, 
				 const double b1, 
				 const double x)
{
	return (1.0 - 1.0 / (1.0 + exp(b0 + b1 * x)));
}

/* 3-param Log-logit curve: b2 / (1 + exp(b0 + b1 * x) ) */
double ls3_fun(const double b0, 
				 const double b1, 
				 const double b2, 
				 const double x)
{
	return (1.0 - b2 / (1.0 + exp(b0 + b1 * x)));
}

void cls3f_error(const size_t n, 
	const double *b, 
	void * fparams, 
	double * fval) 
{
	/*
	Error function of the three parameter dose-response curve 

                    [  /         b2               \ 2 ]
        error = sum [ | -------------------- - y_i |  ]
			     i  [  \1 + exp(b0 + b1*x_i)      /   ]
	*/

	struct ll_param *par = (struct ll_param *) fparams;

	const double * probs = par->probs;
	const double * conc = par->conc;
	const int probs_size = par->probs_size;

	double curry;

	*fval = 0;
	for (int i = 0; i < probs_size; i++) {
		curry = ls3_fun(b[0], b[1], b[2], conc[i]) - probs[i];
		*fval += curry * curry;
	}

}


int ls2_func_f (const gsl_vector * b, 
				void *params, 
				gsl_vector * f)
{
	struct data *d = (struct data *) params;
	double b0 = gsl_vector_get(b, 0);
	double b1 = gsl_vector_get(b, 1);

	for (int i = 0; i < d->n; ++i)
		{
			double x = d->x[i];
			double yi = d->y[i];
			double y = ls2_fun(b0, b1, x);

			gsl_vector_set(f, i, yi - y);
		}
	return GSL_SUCCESS;
}

int ls3_func_f (const gsl_vector * b, 
				void *params, 
				gsl_vector * f)
{
	struct data *d = (struct data *) params;
	double b0 = gsl_vector_get(b, 0);
	double b1 = gsl_vector_get(b, 1);
	double b2 = gsl_vector_get(b, 2);

	for (int i = 0; i < d->n; ++i)
		{
			double x = d->x[i];
			double yi = d->y[i];
			double y = ls3_fun(b0, b1, b2, x);
			gsl_vector_set(f, i, yi - y);
		}
	return GSL_SUCCESS;
}

void callback(const size_t iter, 
				void *params,
				const gsl_multifit_nlinear_workspace *w)
{
	gsl_vector *f = gsl_multifit_nlinear_residual(w);
	gsl_vector *b = gsl_multifit_nlinear_position(w);
	double avratio = gsl_multifit_nlinear_avratio(w);
	double rcond;
	(void) params; /* not used */
	/* compute reciprocal condition number of J(b) */
	gsl_multifit_nlinear_rcond(&rcond, w);

	fprintf(stderr, "iter %2zu: a = %.4f, b = %.4f, c = %.4f, |a|/|v| = %.4f cond(J) = %8.4f, |f(b)| = %.4f\n",
					iter,
					gsl_vector_get(b, 0),
					gsl_vector_get(b, 1),
					gsl_vector_get(b, 2),
					avratio,
					1.0 / rcond,
					gsl_blas_dnrm2(f));
}

void solve_system(gsl_vector *b, 
				gsl_multifit_nlinear_fdf *fdf,
				gsl_multifit_nlinear_parameters *params)
{
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
	const size_t max_iter = 200;
	const double xtol = 1.0e-8;
	const double gtol = 1.0e-8;
	const double ftol = 1.0e-8;
	const size_t n = fdf->n;
	const size_t p = fdf->p;
	gsl_multifit_nlinear_workspace *work =
		gsl_multifit_nlinear_alloc(T, params, n, p);
	gsl_vector * f = gsl_multifit_nlinear_residual(work);
	gsl_vector * y = gsl_multifit_nlinear_position(work);
	int info;
	double chisq0, chisq, rcond;
	/* initialize solver */
	gsl_multifit_nlinear_init(b, fdf, work);
	/* store initial cost */
	gsl_blas_ddot(f, f, &chisq0);

	gsl_multifit_nlinear_driver(max_iter, 
								xtol, gtol, ftol,
								NULL, NULL, 
								&info, work);
	/* store final cost */
	gsl_blas_ddot(f, f, &chisq);
	/* store cond(J(b)) */
	gsl_multifit_nlinear_rcond(&rcond, work);
	gsl_vector_memcpy(b, y);


	gsl_multifit_nlinear_free(work);

}
