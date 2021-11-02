/*
stats/funclib/ll_mod.h
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

/*
 * Contains the log-likelihood functions and jacobians for curve-fitting.
 * Also contains the error function for least squares fitting. 
 */

struct ll_param {
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

void ll_all_AIC(const double *b, 
  void * fparams, 
  double *fval);


// double ls(const double *b,
// 	const double *probs, 
// 	const double *conc, 
// 	const int probs_size, 
// 	const int num_param) ;

// double ls2(const double *b,const double *probs, 
// 	const double *conc, const int probs_size) ;

// double ls3(const double *b,const double *probs, 
// 	const double *conc, const int probs_size) ;
