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


void ll3_min(double *b, //The minimal value that is found.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare, 
	double *beta);

void ll2_min(double *b, //The minimal value that is found.
	double *probs, 
	double *conc, 
	int probs_size, 
	double minimum, //The function value
	double sigsquare);

void ll3_array_min(
	int probs_size, 
	int n_iters,
	double b[][n_iters], //The minima value that are found. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double minimum, //The function value
	double sigsquare, 
	double *beta);

void ll2_array_min(
	int probs_size, 
	int n_iters,
	double b[][n_iters], //The minima value that are found. Should be of size n_iters x 3
	double probs[][probs_size], //Should be of size n_iters x probs_size
	double *conc, 
	double minimum, //The function value
	double sigsquare);