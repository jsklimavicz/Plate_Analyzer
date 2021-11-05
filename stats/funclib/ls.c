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
	return (1.0 / (1.0 + exp(b0 + b1 * x)));
}

/* 3-param Log-logit curve: b2 / (1 + exp(b0 + b1 * x) ) */
double ls3_fun(const double b0, 
				 const double b1, 
				 const double b2, 
				 const double x)
{
	return (b2 / (1.0 + exp(b0 + b1 * x)));
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

int ls2_func_df (const gsl_vector * b, 
				void *params, 
				gsl_matrix * J)
{
	struct data *d = (struct data *) params;
	double b0 = gsl_vector_get(b, 0);
	double b1 = gsl_vector_get(b, 1);

	for (int i = 0; i < d->n; ++i)
		{
			double x = d->x[i];
			double xi = exp(b0 + b1*x);
			double alpha = 1.0 + xi;
			double Da = - xi / (alpha * alpha);

			gsl_matrix_set(J, i, 0, -Da);
			gsl_matrix_set(J, i, 1, -x*Da);
		}
	return GSL_SUCCESS;
}

int ls3_func_df (const gsl_vector * b, 
				void *params, 
				gsl_matrix * J)
{
	struct data *d = (struct data *) params;
	double b0 = gsl_vector_get(b, 0);
	double b1 = gsl_vector_get(b, 1);
	double b2 = gsl_vector_get(b, 2);

	for (int i = 0; i < d->n; ++i)
		{
			double x = d->x[i];
			double xi = exp(b0 + b1*x);
			double alpha = 1.0 + xi;

			double Da = - b2 * xi / (alpha * alpha);
			double Db = x*Da;
			double Dc = 1.0 / alpha;

			gsl_matrix_set(J, i, 0, -Da);
			gsl_matrix_set(J, i, 1, -Db);
			gsl_matrix_set(J, i, 2, -Dc);

		}
	return GSL_SUCCESS;
}

int ls2_func_fvv (const gsl_vector * b, 
				const gsl_vector * v,
				void *params, 
				gsl_vector * fvv)
{
	struct data *d = (struct data *) params;
	double b0 = gsl_vector_get(b, 0);
	double b1 = gsl_vector_get(b, 1);
	double vb0 = gsl_vector_get(v, 0);
	double vb1 = gsl_vector_get(v, 1);
	size_t i;

	for (i = 0; i < d->n; ++i)
		{
			double x = d->x[i];
			double xi = exp(b0 + b1*x);
			double alpha = 1.0 + xi;
			double gamma = xi / (alpha * alpha);

			double Daa = gamma * (xi - 1.0) / alpha;
			double Dab = x * Daa;
			double Dbb = x * Dab;
			double sum;

			sum = vb0 * vb0 * Daa + 2.0 * vb0 * vb1 * Dab +
						  vb1 * vb1 * Dbb;

			gsl_vector_set(fvv, i, -sum);
		}

	return GSL_SUCCESS;
}

int ls3_func_fvv (const gsl_vector * b, 
				const gsl_vector * v,
				void *params, 
				gsl_vector * fvv)
{
	struct data *d = (struct data *) params;
	double b0 = gsl_vector_get(b, 0);
	double b1 = gsl_vector_get(b, 1);
	double b2 = gsl_vector_get(b, 2);
	double vb0 = gsl_vector_get(v, 0);
	double vb1 = gsl_vector_get(v, 1);
	double vb2 = gsl_vector_get(v, 2);
	size_t i;

	for (i = 0; i < d->n; ++i)
		{
			double x = d->x[i];
			double xi = exp(b0 + b1*x);
			double alpha = 1.0 + xi;
			double gamma = xi / (alpha * alpha);

			double Daa = b2 * gamma * (xi - 1.0) / alpha;
			double Dab = x * Daa;
			double Dac = -gamma;
			double Dbb = x * Dab;
			double Dbc = -x * gamma;
			double sum;

			sum = 		  vb0 * vb0 * Daa +
					2.0 * vb0 * vb1 * Dab +
					2.0 * vb0 * vb2 * Dac + 
						  vb1 * vb1 * Dbb +
					2.0 * vb1 * vb2 * Dbc;

			gsl_vector_set(fvv, i, -sum);
			// if (i == 1){printf("DerDev: b0=%f, b1=%f, b2=%f, x=%f \n", b0, b1, b2, x);}
			// if (i == 1){printf("\tv0=%f, v1=%f, v2=%f\n", vb0, vb1, vb2);}
			// if (i == 1){printf("\tDaa=%f, Dab=%f, Dac=%f, Dbb=%f, Dbc=%f\n", Daa, Dab, Dac, Dbb, Dbc);}
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



	/* iterate until convergence */
	// gsl_multifit_nlinear_driver(max_iter, 
	// 							xtol, gtol, ftol,
	// 							callback, NULL, 
	// 							&info, work);
	gsl_multifit_nlinear_driver(max_iter, 
								xtol, gtol, ftol,
								NULL, NULL, 
								&info, work);
	/* store final cost */
	gsl_blas_ddot(f, f, &chisq);
	/* store cond(J(b)) */
	gsl_multifit_nlinear_rcond(&rcond, work);
	gsl_vector_memcpy(b, y);

	/* print summary */

	// fprintf(stderr, "NITER         = %zu\n", gsl_multifit_nlinear_niter(work));
	// fprintf(stderr, "NFEV          = %zu\n", fdf->nevalf);
	// fprintf(stderr, "NJEV          = %zu\n", fdf->nevaldf);
	// fprintf(stderr, "NAEV          = %zu\n", fdf->nevalfvv);
	// fprintf(stderr, "initial cost  = %.12e\n", chisq0);
	// fprintf(stderr, "final cost    = %.12e\n", chisq);
	// fprintf(stderr, "final x       = ( ");
	// for (int ii = 0; ii < p; ii++){
	// 	fprintf(stderr, "%.12e ", gsl_vector_get(b, ii));
	// }
	// fprintf(stderr, ")\n");
	// fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);

	gsl_multifit_nlinear_free(work);

}
