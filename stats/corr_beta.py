# stats/corr_beta.py
# This file is a component of Plate_Analyzer, which can count Drosophila
# L1 larvae, classify them as alive or dead, and determine dose-repsonse
# information based on live/dead count data. 

# Copyright (C) 2021 James Klimavicz

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import scipy.stats as st
import numpy as np


def corr_matrix(plate_ids, rho = 0.10):
	'''
	Generates a len(plate_ids) x len(plate_ids) block matrix with the number of 
	blocks equal to the numebr of unique plates for a compound. 
	Rho is a correlation coefficient for a Gaussian MVN distribution for generating
	correlated beta variables by means of a copula. 
	'''
	mat = np.eye(len(plate_ids))
	for p in set(plate_ids):
		locations = [j for j, x in enumerate(plate_ids) if x == p]
		for i in locations:
			for j in locations:
				if i > j:
					val = rho/(2**(abs(i-j)-1))
					mat[i,j] = val
					mat[j,i] = val
	return mat

def corr_beta_vars(plate_ids, live_count, dead_count, size = 1, scale = 0.5, null_scale = 0.25, rho = 0.10):
	#TODO: Implement Heldane's prior. 
	if null_scale < np.sqrt(np.finfo(float).eps): scale =  np.sqrt(np.finfo(float).eps)
	#Generates the correlated beta variables
	#sanity checks
	assert len(plate_ids) == len(live_count)
	assert len(dead_count) == len(live_count)
	#get a correlation matrix
	corr_mat = corr_matrix(plate_ids, rho = rho)

	#generate multivariate normal random variables
	norm_var = np.random.default_rng().multivariate_normal(np.array([0.] *len(plate_ids)), corr_mat, size = size)
	#find normal cdf values 
	norm_prob = st.norm.cdf(norm_var)
	#generate beta variables from cdf values
	scale_add = np.where(np.logical_or(dead_count == 0, live_count == 0), null_scale, scale)
	beta_var = st.beta.ppf(norm_prob, dead_count + scale_add, live_count + scale_add)
	return beta_var
