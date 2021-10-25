# stats/functionfit.py
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

from ctypes import *
import numpy as np
import os
import platform
import stats.utils as utils

class FunctionFit():
	SS = 1.0e6
	BP=np.array([1.5,1.01])

	def __init__(self):
		p = platform.platform()
		if "linux" in p.lower():
			lib_path = "./stats/funclib/cloglik.so"
		elif "windows" in p.lower():
			lib_path = "./stats/funclib/cloglik.dll"
		if os.path.exists(lib_path):
			self.cloglik = CDLL(lib_path)
			self.use_C_lib = True
			#Set LL3 vars
			self.ll3c = self.cloglik.ll3_min
			self.ll3c.argtypes = (np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
					np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
					np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
					c_int, 
					c_double, 
					c_double, 
					np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'))
			#Set LL2 vars
			self.ll2c = self.cloglik.ll2_min
			self.ll2c.argtypes = (np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
					np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
					np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
					c_int, 
					c_double, 
					c_double)
		else:
			self.cloglik = None
			self.use_C_lib = False

	def min_ll3(self, b, probs, conc):
		if self.use_C_lib:
			funmin = 0
			self.ll3c(b, 
					probs, 
					conc, 
					len(probs), 
					funmin,
					self.SS,
					self.BP)
			return b
		else:
			res = minimize(self.__ll3p, b, args = (probs, self.conc), 
					method = 'BFGS', jac = __ll3p_jac)
			return res.x

	def min_ll2(self, b, probs, conc):
		if self.use_C_lib:
			funmin = 0
			self.ll2c(b, 
					probs, 
					conc, 
					len(probs), 
					funmin,
					self.SS)
			return b
		else:
			res = minimize(self.__ll2p, b, args = (probs, self.conc), 
					method = 'BFGS', jac = __ll2p_jac)
			return res.x

	@staticmethod
	@utils.surpress_warnings
	def __ll3p(b, probs, conc, sigma_squared = 1e6, beta_param=[1.5,1.01]):
		'''
		Log-likelihood function of the three parameter dose-response curve 
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
		b2 ~ Beta(beta_param).
		'''
		b0, b1, b2 = b
		if (b2 <= 1e-10 ): return(1e10)
		xi = np.exp(b0+b1*conc)
		alpha = 1+xi # >1.0
		# l = np.log(alpha)
		if (min(1+xi)-b2 <= 1e-10): return(1e10)
		ba = beta_param[0]
		bb = beta_param[1]
		#note: from benchmarking, b0*b0 is approximately 3x faster than b0**2 for a float.
		#MVN Prior
		ll = -(b0*b0 + b1*b1)/(2*sigma_squared) 

		#Beta Prior
		ll += (ba-1)*math.log(b2) + (bb-1)*math.log(1-b2)

		#terms
		ll += sum(probs*np.log(alpha-b2) - np.log(alpha) + (1-probs)*math.log(b2))
		return(-ll)

	@staticmethod
	@utils.surpress_warnings
	def __ll3p_jac(b, probs, conc,sigma_squared = 1e6, weibull_param=[2,1]):
		'''
		Jacobian of the log-likelihood function of the three 
		parameter dose-response curve 
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
		b2 ~ Beta(beta_param).
		'''
		b0, b1, b2 = b

		ba = beta_param[0]
		bb = beta_param[1]

		xi = np.exp(b0+b1*conc)
		alpha = 1+xi
		d = (alpha - b2)

		g0 = -b0/sigma_squared + sum(probs * xi/( d) - xi/(alpha))
		g1 = -b1/sigma_squared + sum(conc * xi *(probs/(d) - 1/(alpha)) )
		g2 = (ba-1)/b2 - (bb-1)/(1-b2) + sum(-probs/(d)) + sum((1-probs)/b2)
		return(np.array([-g0,-g1,-g2]))


	@staticmethod
	@utils.surpress_warnings
	def __ll2p(b, probs, conc, sigma_squared = 1e6):
		'''
		Log-likelihood function of the two parameter dose-response curve 
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein prior is b0, b1 ~ MVN(0, sigma*I2).
		'''
		b0, b1 = b
		p_sum = sum(probs)
		p_conc_sum = sum(probs*conc)
		ll = -(b0*b0 + b1*b1)/(2*sigma_squared) + b0*p_sum + b1*p_conc_sum - \
				sum(np.log(1 + np.exp(b0+b1*conc)))
		return(-ll)

	@staticmethod
	@utils.surpress_warnings
	def __ll2p_jac(b, probs, conc, sigma_squared = 1e6):
		'''
		Jacobian of the log-likelihood function of the 
		two parameter dose-response curve
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein prior is b0, b1 ~ MVN(0, sigma*I2).
		'''
		b0, b1 = b
		xi = np.exp(b0+b1*conc)
		l = xi/(xi+1)
		g1 = -b0/sigma_squared + sum(probs) - sum(l)
		g2 = -b1/sigma_squared + sum(conc*probs) - sum(conc*l)
		return(np.array([-g1,-g2]))