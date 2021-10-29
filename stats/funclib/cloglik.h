/*
cloglik.h
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


void ll3_min(double *b, //The val to be minimized.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare, 
	double *beta);

void ll2_min(double *b, //The val to be minimized.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare);

void ll3_array_min(
	int probs_size, 
	int n_iters,
	double b[][3], //The vals to be minimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double *minimum, //The function values
	double sigsquare, 
	double *beta);

void ll2_array_min(
	int probs_size, 
	int n_iters,
	double b[][2], //The vals to be minimized. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double *minimum, //The function values
	double sigsquare);

/*
This function for for calculating both the ll2 curve and the ll3 curve, and
then the AIC is used to determine whether the 2-parameter or 3-parameter curve
should be used. 
*/
void ll2_ll3_AIC(double *b, // The vals to be minimized. Must have length 3.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare,
	double *beta);

void array_ll2_ll3_AIC(
	int probs_size, 
	int n_iters,
	double b[][3], //The minima value that are found. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double *minimum, //The function value
	double sigsquare, 
	double *beta);

void ls_driver(
	const int nparam, //number of parameters to use
	const int n, //number of data points
	double *user_b, //initial guess for b (and output). Must be len 3
	double *x, //concentration values
	double *y,// y-values
	const int method); 

void ls_array_min(
	const int nparam,
	const int probs_size, 
	const int n_iters,
	double b[][3], //The vals to be minimized. Should be of size n_iters x 3
	double *conc,
	double probs[][probs_size], //Should be of size n_iters x probs_size
	const int method);
