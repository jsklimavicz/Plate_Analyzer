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
import ctypes
import numpy as np
import os
import platform
import stats.utils as utils
import math
from scipy.optimize import minimize, least_squares

class FunctionFit():
	SS = 1.0e6
	BP=np.array([1.5,1.01])
	'''
	Least squares methods:
	LS_METHOD = 0: Levenberg-Marquardt
	LS_METHOD = 1: Levenberg-Marquardt with geodesic acceleration
	LS_METHOD = 2: Dogleg
	LS_METHOD = 3: Double dogleg
	LS_METHOD = 4: 2D subspace
	'''
	LS_METHOD = 4
	OPTIM_METHOD = 5

	def __init__(self, **kwargs):
		self.set_default_params(**kwargs)
		p = platform.platform()
		if "linux" in p.lower():
			lib_path = "./stats/funclib/cloglik.so"
			func = CDLL
		elif "windows" in p.lower():
			lib_path = "./stats/funclib/cloglik.dll"
			func = ctypes.WinDLL
			#USER MUST UPDATE THIS WITH THE APPROPRIATE PATH TO THE DIRECTORY
			#CONTAINING THE libgsl.dll FILE AND OTHER APPROPRIATE .dll FILES.
			os.add_dll_directory("C:/msys64/mingw64/bin")
		else:
			lib_path = None
		if os.path.exists(lib_path) and 'USE_CLIB' in kwargs and kwargs.get('USE_CLIB'):

			#necessary for the case when the library loads but is not working,
			#e.g. because gsl is on not the system. 
			try:
				self.cloglik = func(lib_path)
				#Single shot optimizers
				single_min_params = [np.ctypeslib.ndpointer(dtype=np.float64,
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, 
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, 
												ndim=1, flags='C_CONTIGUOUS'),
						c_int, 
						c_double, 
						c_double, 
						np.ctypeslib.ndpointer(dtype=np.float64, 
												ndim=1, flags='C_CONTIGUOUS'),
						c_int] #OPTIM_METHOD
				#LL3 minimzer
				self.ll3c = self.cloglik.ll3_min
				self.ll3c.argtypes = (single_min_params)

				#LL2 minimzer
				self.ll2c = self.cloglik.ll2_min
				#all but beta params
				self.ll2c.argtypes = (*single_min_params[:-2],single_min_params[-1])

				#AIC-based minimizer
				self.ll23cAIC = self.cloglik.ll2_ll3_AIC
				self.ll23cAIC.argtypes = (single_min_params)

				#Array-based optimizers
				array_min_params = [c_int, #number of probs/trial
						c_int, #number of iters
						np.ctypeslib.ndpointer(dtype=np.float64, #minima
												ndim=2, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #prob array
												ndim=2, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #conc
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #func vals
												ndim=1, flags='C_CONTIGUOUS'),
						c_double, #sigma**2
						np.ctypeslib.ndpointer(dtype=np.float64, #beta perams
												ndim=1, flags='C_CONTIGUOUS'),
						c_int] #OPTIM_METHOD
				#LL3 Array minimizer
				self.ll3ca = self.cloglik.ll3_array_min
				self.ll3ca.argtypes = (array_min_params)

				#LL2 Array minimizer
				self.ll2ca = self.cloglik.ll2_array_min
				#all but beta params
				self.ll2ca.argtypes = (*array_min_params[:-2],array_min_params[-1]) 

				#AIC-based best LL picker Array minimizer
				self.ll23aAIC = self.cloglik.array_ll2_ll3_AIC
				self.ll23aAIC.argtypes = (array_min_params)

				#Least-Squares fitters
				#Single Curve
				self.lsfit = self.cloglik.ls_driver
				self.lsfit.argtypes = (c_int, #number of parameters (2 or 3)
				c_int, # number of data points
				np.ctypeslib.ndpointer(dtype=np.float64, #input b
											ndim=1, flags='C_CONTIGUOUS'),
				np.ctypeslib.ndpointer(dtype=np.float64, #conc values
											ndim=1, flags='C_CONTIGUOUS'),
				np.ctypeslib.ndpointer(dtype=np.float64, #prob values
											ndim=1, flags='C_CONTIGUOUS'))

				#Array
				self.lsarray = self.cloglik.ls_array_min
				self.lsarray.argtypes = (c_int, #number of parameters (2 or 3)
				c_int, # number of data points
				c_int, # n_iters
				np.ctypeslib.ndpointer(dtype=np.float64, #input b
											ndim=2, flags='C_CONTIGUOUS'),
				np.ctypeslib.ndpointer(dtype=np.float64, #conc values
											ndim=1, flags='C_CONTIGUOUS'),
				np.ctypeslib.ndpointer(dtype=np.float64, #prob values
											ndim=2, flags='C_CONTIGUOUS'))

				self.use_C_lib = True
			except:
				self.use_C_lib =  False
				self.cloglik = None
		else:
			self.cloglik = None
			self.use_C_lib = False



	def set_default_params(self, **kwargs):
		self.use_C_lib = kwargs.get('USE_CLIB') if 'USE_CLIB' in kwargs else False
		self.LS_METHOD = kwargs.get('LS_METHOD') if 'LS_METHOD' in kwargs else 4
		self.OPTIM_METHOD = kwargs.get('OPTIM_METHOD') if 'OPTIM_METHOD' in kwargs else 2
		self.SS = kwargs.get('LL_SIGMA')**2 if 'LL_SIGMA' in kwargs else 1e6
		self.BP=np.array([1.5,1.01])
		self.BP[0] = kwargs.get('LL_BETA1') if 'LL_BETA1' in kwargs else 1.5
		self.BP[1] = kwargs.get('LL_BETA2') if 'LL_BETA2' in kwargs else 1.01

	def array_ll3(self, b_input, prob_array, conc_list, 
							sigma_squared = None, beta_param=None,
							optim_method = None):
		if not sigma_squared: sigma_squared = self.SS
		if not beta_param: beta_param = self.BP
		if not optim_method: optim_method = self.OPTIM_METHOD

		niters, prob_ct = prob_array.shape
		if self.use_C_lib:
			funmin = np.zeros(niters)
			self.ll3ca(prob_ct,
						niters,
						b_input, #must be size niters * 3
						prob_array, #must be size niters * prob_ct
						conc_list, #must be length prob_ct
						funmin,
						sigma_squared,
						beta_param,
						optim_method)
		else:
			for i in range(niters):
				b_input[i] = self.min_ll3(b_input[i], prob_array[i],
									conc_list, sigma_squared, beta_param)
		return b_input

	def array_ll2(self, b_input, prob_array, conc_list, sigma_squared = None,
							optim_method = None):

		if not sigma_squared: sigma_squared = self.SS
		if not optim_method: optim_method = self.OPTIM_METHOD

		# print(b_input)
		print(b_input.shape)
		# print(prob_array)
		print(prob_array.shape)
		# print(conc_list)
		print(conc_list.shape)
		print(sigma_squared)
		print(optim_method)
		print(self.use_C_lib)

		niters, prob_ct = prob_array.shape
		if self.use_C_lib: 
			funmin = np.zeros(niters)
			self.ll2ca(prob_ct,
						niters,
						b_input, #must be size niters * 2
						prob_array, #must be size niters * prob_ct
						conc_list, #must be length prob_ct
						funmin,
						sigma_squared,
						optim_method)
		else:
			b_out = np.zeros_like(b_input)
			for i in range(niters):
				b_input[i] = self.min_ll2(b_input[i], prob_array[i], 
									conc_list, sigma_squared)
		return b_input

	def array_ll23AIC(self, b_input, prob_array, conc_list, 
							sigma_squared = None, beta_param=None,
							optim_method = None):
		if not sigma_squared: sigma_squared = self.SS
		if not beta_param: beta_param = self.BP
		if not optim_method: optim_method = self.OPTIM_METHOD
		niters, prob_ct = prob_array.shape
		if self.use_C_lib: 
			funmin = np.zeros(niters)
			self.ll23aAIC(prob_ct,
						niters,
						b_input, #must be size niters * 3
						prob_array, #must be size niters * prob_ct
						conc_list, #must be length prob_ct
						funmin,
						sigma_squared,
						beta_param,
						optim_method)
		else:
			for i in range(niters):
				b_input[i] = self.min_llAIC(b_input[i], prob_array[i], 
									conc_list, sigma_squared, beta_param)
		return b_input

	def array_ls(self, nparam, b_input, prob_array, conc_list, ls_method = None):
		if not ls_method: ls_method = self.LS_METHOD
		if self.use_C_lib:
			niters, prob_ct = prob_array.shape
			self.lsarray(nparam, 
					prob_ct, 
					niters,
					b_input, 
					conc_list, 
					prob_array, 
					ls_method)
		else:
			for i in range(niters):
				b_input[i] = FunctionFit.ll_ls(b_input[i], nparam, prob_array[i], conc_list)
		return b_input

	def min_ls(self, b, nparam, probs, conc, ls_method = None):
		if not ls_method: ls_method = self.LS_METHOD
		if self.use_C_lib:
			nprob = len(probs)
			self.lsfit(nparam, 
					nprob, 
					b, 
					conc, 
					probs, 
					ls_method)
			return b
		else:
			return FunctionFit.ll_ls(b, nparam, probs, conc)

	def min_ll3(self, b, probs, conc, sigma_squared = None, beta_param=None,
					optim_method = None):
		if not sigma_squared: sigma_squared = self.SS
		if not beta_param: beta_param = self.BP
		if not optim_method: optim_method = self.OPTIM_METHOD
		if self.use_C_lib:
			funmin = 0
			self.ll3c(b, 
					probs, 
					conc, 
					len(probs), 
					funmin,
					sigma_squared,
					beta_param,
					optim_method)
			return b
		else:
			res = minimize(FunctionFit.ll3p, b, 
					args = (probs, conc, sigma_squared, beta_param), 
					method = 'BFGS', jac = FunctionFit.ll3p_jac)
			return res.x

	def min_llAIC(self, b, probs, conc, sigma_squared = None, beta_param=None,
					optim_method = None):
		if not sigma_squared: sigma_squared = self.SS
		if not beta_param: beta_param = self.BP
		if not optim_method: optim_method = self.OPTIM_METHOD
		if self.use_C_lib:
			funmin = 0
			self.ll23cAIC(b, 
					probs, 
					conc, 
					len(probs), 
					funmin,
					sigma_squared,
					beta_param,
					optim_method)
			return b
		else:
			return self.ll23AIC_min(b, probs, conc, sigma_squared, beta_param)

	def min_ll2(self, b, probs, conc, sigma_squared = None, optim_method = None):
		if not sigma_squared: sigma_squared = self.SS
		if not optim_method: optim_method = self.OPTIM_METHOD

		if self.use_C_lib:
			funmin = 0
			self.ll2c(b, 
					probs, 
					conc, 
					len(probs), 
					funmin,
					sigma_squared)
			return b
		else:
			res = minimize(FunctionFit.ll2p, b, 
					args = (probs, conc, sigma_squared), 
					method = 'BFGS', jac = FunctionFit.ll2p_jac)
			return res.x

	@staticmethod
	@utils.surpress_warnings
	def ll_ls(b, nparam, probs, conc):
		return least_squares(FunctionFit.least_squares_error,
							b, args=(nparam, probs, conc)).x

	@staticmethod
	@utils.surpress_warnings
	def ll3p(b, probs, conc, sigma_squared = SS, beta_param=BP):
		'''
		Log-likelihood function of the three parameter dose-response curve 
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
		b2 ~ Beta(beta_param).
		'''
		b0, b1, b2 = b
		if ((b2 <= 1e-10) or (b2 >=1) ): return(1e10)
		xi = np.exp(b0+b1*conc)
		alpha = 1+xi # >1.0
		# l = np.log(alpha)
		if (min(1+xi)-b2 <= 1e-10): return(1e10)
		ba = beta_param[0]
		bb = beta_param[1]
		#MVN Prior
		ll = -(b0*b0 + b1*b1)/(2*sigma_squared) 
		#Beta Prior
		ll += (ba-1)*math.log(b2) + (bb-1)*math.log(1-b2)
		#terms
		ll += sum(probs*np.log(alpha-b2) - np.log(alpha) + \
					(1-probs)*math.log(b2))
		return(-ll)

	@staticmethod
	@utils.surpress_warnings
	def ll_for_AIC(b, probs, conc, sigma_squared = SS, beta_param=BP):
		b0, b1, b2 = b

		#MVN Prior
		ll = -(b0*b0 + b1*b1)/(2*sigma_squared) 

		if (abs(b2 - 1) < 1e-6): #ll2
			p_sum = sum(probs)
			p_conc_sum = sum(probs*conc)
			ll = -(b0*b0 + b1*b1)/(2*sigma_squared) + b0*p_sum + b1*p_conc_sum - \
					sum(np.log(1 + np.exp(b0+b1*conc)))
		else: #ll3
			xi = np.exp(b0+b1*conc)
			alpha = 1+xi
			ll += (beta_param[0]-1)*math.log(b2) + (beta_param[1]-1)*math.log(1-b2) #Beta Prior
			ll += sum(probs*np.log(alpha-b2) - np.log(alpha) + \
					(1-probs)*math.log(b2))
		return ll

	@staticmethod
	@utils.surpress_warnings
	def ll23AIC_min(b, probs, conc, sigma_squared = SS, beta_param=BP):
		#fit ll2
		res2 = minimize(FunctionFit.ll2p, b[0:2], 
					args = (probs, conc, sigma_squared), 
					method = 'BFGS', jac = FunctionFit.ll2p_jac)
		b2 = np.array([*res2.x, 1.0])
		ll2 = FunctionFit.ll_for_AIC(b2, probs, conc, sigma_squared)
		#fit ll3
		res3 = minimize(FunctionFit.ll3p, b, 
					args = (probs, conc, sigma_squared, beta_param), 
					method = 'BFGS', jac = FunctionFit.ll3p_jac)
		ll3 = FunctionFit.ll_for_AIC(res3.x, probs, conc, sigma_squared, beta_param)
		AIC2 = 4 - 2*ll2
		AIC3 = 6 - 2*ll3

		if (AIC2 < AIC3):
			return b2
		else:
			return res3.x

	@staticmethod
	@utils.surpress_warnings
	def ll3p_jac(b, probs, conc,sigma_squared = SS, beta_param=BP):
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
	def ll2p(b, probs, conc, sigma_squared = SS):
		'''
		Log-likelihood function of the two parameter dose-response curve 
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein prior is b0, b1 ~ MVN(0, sigma*I2).
		'''
		b0 = b[0]
		b1 = b[1]
		p_sum = sum(probs)
		p_conc_sum = sum(probs*conc)
		ll = -(b0*b0 + b1*b1)/(2*sigma_squared) + b0*p_sum + b1*p_conc_sum - \
				sum(np.log(1 + np.exp(b0+b1*conc)))
		return(-ll)

	@staticmethod
	@utils.surpress_warnings
	def ll2p_jac(b, probs, conc, sigma_squared = SS):
		'''
		Jacobian of the log-likelihood function of the 
		two parameter dose-response curve
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein prior is b0, b1 ~ MVN(0, sigma*I2).
		'''
		b0 = b[0]
		b1 = b[1]
		xi = np.exp(b0+b1*conc)
		l = xi/(xi+1)
		g1 = -b0/sigma_squared + sum(probs) - sum(l)
		g2 = -b1/sigma_squared + sum(conc*probs) - sum(conc*l)
		return(np.array([-g1,-g2]))

	@staticmethod
	@utils.surpress_warnings
	def least_squares_error(b, nparam, probs, conc):
		'''
		Dose-response curve for least-squares fitting. If len(b)=2, then 
		b2 = 1; else if len(b)=3, then b2 = b[2]
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		'''
		if nparam == 2:
			return np.array((1-1/(1 + np.exp(b[0] + b[1]*conc)) - probs)**2)
		else:
			if b[2] > 1:
				return 1e10
			else:
				return np.array((1-b[2]/(1 + np.exp(b[0] + b[1]*conc)) - probs)**2)