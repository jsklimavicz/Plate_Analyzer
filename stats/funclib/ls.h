/*
stats/funclib/ls.h
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

void cls3f_error(const size_t n, 
        const double *b, 
        void * fparams, 
        double * fval);

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
