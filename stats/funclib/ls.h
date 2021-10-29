#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>

struct data{
  double *x;
  double *y;
  size_t n;
};

double ls2_fun(const double b0, 
         const double b1, 
         const double x);

double ls3_fun(const double b0, 
         const double b1, 
         const double b2, 
         const double x);

int ls2_func_f (const gsl_vector * b, 
        void *params, 
        gsl_vector * f);

int ls3_func_f (const gsl_vector * b, 
        void *params, 
        gsl_vector * f);

int ls2_func_df (const gsl_vector * b, 
        void *params, 
        gsl_matrix * J);

int ls3_func_df (const gsl_vector * b, 
        void *params, 
        gsl_matrix * J);

int ls2_func_fvv (const gsl_vector * b, 
        const gsl_vector * v,
        void *params, 
        gsl_vector * fvv);

int ls3_func_fvv (const gsl_vector * b, 
        const gsl_vector * v,
        void *params, 
        gsl_vector * fvv);

void solve_system(gsl_vector *b, 
        gsl_multifit_nlinear_fdf *fdf,
        gsl_multifit_nlinear_parameters *params);