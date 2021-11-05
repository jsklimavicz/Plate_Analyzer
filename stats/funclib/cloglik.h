/*
stats/funclib/cloglik.h
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

//SINGLE-FIT METHODS//

void ll2_min(double *b, //Vector of parameters for optimization
	double *probs, //Array of probs
	double *conc, //Array of concentrations
	int probs_size, //Length of the prob/conc arrays
	double minimum, //The function value
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	const int method); //ID of method to be usedfor optimization.

void ll3_min(double *b, //Vector of parameters for optimization
	double *probs, //Array of probs
	double *conc, //Array of concentrations
	int probs_size, //Length of the prob/conc arrays
	double minimum, //The function value
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	double *beta, //2-value array with parameters for beta prior on b2
	const int method); //ID of method to be usedfor optimization. 

void ll2_ll3_AIC(double *b, //Vector of parameters for optimization
	double *probs, //Array of probs
	double *conc, //Array of concentrations
	int probs_size, //Length of the prob/conc arrays
	double minimum, //The function value
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	double *beta, //2-value array with parameters for beta prior on b2
	const int method); //ID of method to be usedfor optimization. 

void ls_driver(
	const int nparam, //number of parameters to use in dose-response curve
	const int n, //number of data points
	double *user_b, //initial guess for b (and output). Must be len 3
	double *x, //concentration values
	double *y, //probability values
	const int method); //least squares method to use

//ARRAY-BASED METHODS//

void ll2_array_min(
	int probs_size, //Number of vals in the prob/conc arrays
	int n_iters, //Number of different sets of probs
	double b[][2], //The vals to be optimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, //Array of concentrations
	double *minimum, //The function values
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	const int method); //ID of method to be usedfor optimization.

void ll3_array_min(
	int probs_size, //Number of vals in the prob/conc arrays
	int n_iters, //Number of different sets of probs
	double b[][3], //The vals to be optimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, //Array of concentrations
	double *minimum, //The function values
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	double *beta, //2-value array with parameters for beta prior on b2
	const int method); //ID of method to be usedfor optimization. 

void array_ll2_ll3_AIC(
	int probs_size, //Number of vals in the prob/conc arrays
	int n_iters, //Number of different sets of probs
	double b[][3], //The vals to be optimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, //Array of concentrations
	double *minimum, //The function values
	double *lb, //Lower bound for constrained optim
	double *ub, //Upper bound for constrained optim
	double sigsquare, //Value of sigma**2 for MVN prior on b0 and b1
	double *beta, //2-value array with parameters for beta prior on b2
	const int method); //ID of method to be usedfor optimization. 

void ls_array_min(
	const int nparam, //number of parameters to use in dose-response curve
	const int probs_size, //Number of vals in the prob/conc arrays
	const int n_iters, //Number of different sets of probs
	double b[][3], //The vals to be optimized. Should be of size n_iters x 3
	double *conc, //concentration values
	double probs[][probs_size], //Should be of size n_iters x probs_size
	const int method); //least squares method to use

void cls3_min(double *b, //Vector of parameters for optimization
	double *probs, //Array of probs
	double *conc, //Array of concentrations
	int probs_size, //Length of the prob/conc arrays
	double minimum, //The function value
	double *lb, //Lower bound for constrained optim
	double *ub); //Upper bound for constrained optim 

void cls3_array_min(
	int probs_size, //Number of vals in the prob/conc arrays
	int n_iters, //Number of different sets of probs
	double b[][3], //The vals to be optimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, //Array of concentrations
	double *minimum, //The function values
	double *lb, //Lower bound for constrained optim
	double *ub); //Upper bound for constrained optim 