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
	LB = np.array([-10, 0.05, 0])
	UB = np.array([10, 5, 1])


	def __init__(self, **kwargs):
		self.set_default_params(**kwargs)
		#First check to see if a C shared-object or dynamic linked library
		#exists, and if it does, try using the C library. 
		p = platform.platform()
		if ("linux" in p.lower()) or ("macos" in p.lower()):
			lib_path = "./stats/funclib/cloglik.so"
			func = CDLL
		elif "windows" in p.lower():
			lib_path = "./stats/funclib/cloglik.dll"
			func = ctypes.WinDLL
			#USER MUST UPDATE THIS WITH THE APPROPRIATE PATH TO THE DIRECTORY
			#CONTAINING THE libgsl.dll FILE AND OTHER APPROPRIATE .dll FILES.
			os.add_dll_directory("C:/msys64/mingw64/bin")
		else:
			lib_path = ""
		if os.path.exists(lib_path) and 'USE_CLIB' in kwargs and kwargs.get('USE_CLIB'):
			'''	
			try is necessary for the case when the library loads but is not working,
			e.g. because gsl is on not the system, or it can;t be found, or because 
			the libray was not properly compiled. 
			'''
			try:
				self.__cloglik = func(lib_path)
				#Single shot optimizers
				single_min_params = [
						np.ctypeslib.ndpointer(dtype=np.float64, #vector param
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #probs
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #conc
												ndim=1, flags='C_CONTIGUOUS'),
						c_int, #nparam
						c_double, #funtion min
						np.ctypeslib.ndpointer(dtype=np.float64, #lb
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #ub
												ndim=1, flags='C_CONTIGUOUS'),
						c_double, #sigma_squared
						np.ctypeslib.ndpointer(dtype=np.float64, #beta_param
												ndim=1, flags='C_CONTIGUOUS'),
						c_int] #OPTIM_METHOD
				#LL3 minimzer
				self.__ll3c = self.__cloglik.ll3_min
				self.__ll3c.argtypes = (single_min_params)

				#LL2 minimzer
				self.__ll2c = self.__cloglik.ll2_min
				#all but beta params
				self.__ll2c.argtypes = (*single_min_params[:-2],single_min_params[-1])

				#AIC-based minimizer
				self.__ll23cAIC = self.__cloglik.ll2_ll3_AIC
				self.__ll23cAIC.argtypes = (single_min_params)

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
						np.ctypeslib.ndpointer(dtype=np.float64, #lb
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #ub
												ndim=1, flags='C_CONTIGUOUS'),
						c_double, #sigma**2
						np.ctypeslib.ndpointer(dtype=np.float64, #beta perams
												ndim=1, flags='C_CONTIGUOUS'),
						c_int] #OPTIM_METHOD
				#LL3 Array minimizer
				self.__ll3ca = self.__cloglik.ll3_array_min
				self.__ll3ca.argtypes = (array_min_params)

				#LL2 Array minimizer
				self.__ll2ca = self.__cloglik.ll2_array_min
				#all but beta params
				self.__ll2ca.argtypes = (*array_min_params[:-2],array_min_params[-1]) 

				#AIC-based best LL picker Array minimizer
				self.__ll23aAIC = self.__cloglik.array_ll2_ll3_AIC
				self.__ll23aAIC.argtypes = (array_min_params)

				#Least-Squares fitters
				#Single Curve
				self.__lsfit = self.__cloglik.ls_driver
				self.__lsfit.argtypes = (c_int, #number of parameters (2 or 3)
				c_int, # number of data points
				np.ctypeslib.ndpointer(dtype=np.float64, #input b
											ndim=1, flags='C_CONTIGUOUS'),
				np.ctypeslib.ndpointer(dtype=np.float64, #conc values
											ndim=1, flags='C_CONTIGUOUS'),
				np.ctypeslib.ndpointer(dtype=np.float64, #prob values
											ndim=1, flags='C_CONTIGUOUS'))

				#Array
				self.__lsarray = self.__cloglik.ls_array_min
				self.__lsarray.argtypes = (c_int, #number of parameters (2 or 3)
				c_int, # number of data points
				c_int, # n_iters
				np.ctypeslib.ndpointer(dtype=np.float64, #input b
											ndim=2, flags='C_CONTIGUOUS'),
				np.ctypeslib.ndpointer(dtype=np.float64, #conc values
											ndim=1, flags='C_CONTIGUOUS'),
				np.ctypeslib.ndpointer(dtype=np.float64, #prob values
											ndim=2, flags='C_CONTIGUOUS'))


				#Constrained Least Squares
				self.__cls = self.__cloglik.cls3_min
				self.__cls.argtypes = [
						np.ctypeslib.ndpointer(dtype=np.float64, #minimum
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #prob array
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #conc
												ndim=1, flags='C_CONTIGUOUS'),
						c_int, #number of probabilities
						c_double, #function value
						np.ctypeslib.ndpointer(dtype=np.float64, #lb
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #ub
												ndim=1, flags='C_CONTIGUOUS')]


				#Constrained Least Squares by array
				self.__cls_array = self.__cloglik.cls3_array_min
				self.__cls_array.argtypes = [c_int, #number of probs/trial
						c_int, #number of iters
						np.ctypeslib.ndpointer(dtype=np.float64, #minima
												ndim=2, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #prob array
												ndim=2, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #conc
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #func vals
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #lb
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #ub
												ndim=1, flags='C_CONTIGUOUS')]

				self.use_C_lib = True
			except:
				self.use_C_lib =  False
				self.__cloglik = None
		else:
			self.__cloglik = None
			self.use_C_lib = False
		

	def check_Clib_running(self):
		return self.use_C_lib

	def set_default_params(self, **kwargs):
		self.use_C_lib = kwargs.get('USE_CLIB') if 'USE_CLIB' in kwargs else False
		self.LS_METHOD = kwargs.get('LS_METHOD') if 'LS_METHOD' in kwargs else 4
		self.OPTIM_METHOD = kwargs.get('OPTIM_METHOD') if 'OPTIM_METHOD' in kwargs else 2
		self.SS = kwargs.get('LL_SIGMA')**2 if 'LL_SIGMA' in kwargs else 1e6
		self.BP=np.array([1.5,1.01])
		self.BP[0] = kwargs.get('LL_BETA1') if 'LL_BETA1' in kwargs else 1.5
		self.BP[1] = kwargs.get('LL_BETA2') if 'LL_BETA2' in kwargs else 1.01
		self.LB = kwargs.get('LB') if 'LB' in kwargs else np.array([-10, 0.05, 0])
		self.UB = kwargs.get('UB') if 'UB' in kwargs else np.array([10, 5, 1])
 
	def array_ll2(self, b_input, prob_array, conc_list, sigma_squared = None,
							optim_method = None, lb = None, ub = None):
		'''
		Used to pass an array of probabilities to the C-based or python-based 
			ll2 solver. For data that has prob_ct different concentration/
			probability pairs and niter different iterations of bootstrapped 
			probabilities:
		b_input: 2D np.array of size niters * 2. Should initially contain the
			guesses for b, and will ultimately contain the optimized values.
		prob_array: 2D np.array of size niters * prob_ct containing the 
			bootstrapped probability values.
		conc_list: 1D np.array of length prob_ct. Contains the concentrations
			that are matched to the probabilities. 
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b[0] and b[1]. 
		optim_method: An int from {0,1,2,3,4,5,6,7} to choose the solver. See
			the README for more details. 

		returns: 2D np.array of size niters * 2 containing the optimized 
			values for each bootstrap. 
		'''
		if sigma_squared is None: sigma_squared = self.SS
		if optim_method is None: optim_method = self.OPTIM_METHOD
		if lb is None: lb = self.LB
		if ub is None: ub = self.UB

		niters, prob_ct = prob_array.shape
		if self.use_C_lib: # pass off to the C function if it exists.
			funmin = np.zeros(niters)
			self.__ll2ca(prob_ct, niters, b_input, prob_array, conc_list, 
						funmin, lb[0:2], ub[0:2], sigma_squared, optim_method)
		else:
			for i in range(niters):
				b_input[i] = self.min_ll2(b_input[i], prob_array[i], 
									conc_list, sigma_squared)
		#need to make an niters * 3 array
		return np.c_[b_input,np.ones(niters)]

	def array_ll3(self, b_input, prob_array, conc_list, 
							sigma_squared = None, beta_param=None,
							optim_method = None, lb = None, ub = None):
		'''
		Used to pass an array of probabilities to the C-based or python-based
			ll3 solver. For data that has prob_ct different concentration/
			probability pairs and niter different iterations of bootstrapped 
			probabilities:
		b_input: 2D np.array of size niters * 3. Should initially contain the
			guesses for b, and will ultimately contain the optimized values.
		prob_array: 2D np.array of size niters * prob_ct containing the 
			bootstrapped probability values.
		conc_list: 1D np.array of length prob_ct. Contains the concentrations
			that are matched to the probabilities. 
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b[0] and b[1]. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b[2].
		optim_method: An int from {0,1,2,3,4,5,6,7} to choose the solver. See
			the README for more details. 

		returns: 2D np.array of size niters * 3 containing the optimized 
			values for each bootstrap. 
		'''
		if sigma_squared is None: sigma_squared = self.SS
		if beta_param is None: beta_param = self.BP
		if optim_method is None: optim_method = self.OPTIM_METHOD
		if lb is None: lb = self.LB
		if ub is None: ub = self.UB

		niters, prob_ct = prob_array.shape
		if self.use_C_lib: # pass off to the C function if it exists.
			#funmin is for storing the function minimum, and can be returned
			#to provide the minimum value if needed. 
			funmin = np.zeros(niters)
			self.__ll3ca(prob_ct, niters, b_input, prob_array, conc_list, 
						funmin, lb, ub, sigma_squared, 
						beta_param, optim_method)
		else:
			for i in range(niters):
				b_input[i] = self.min_ll3(b_input[i], prob_array[i],
									conc_list, sigma_squared, beta_param)
		return b_input

	def array_ll23AIC(self, b_input, prob_array, conc_list, 
							sigma_squared = None, beta_param=None,
							optim_method = None, lb = None, ub = None):
		'''
		Used to pass an array of probabilities to the C-based or python-based 
			ll2/ll3 solver that then uses the AIC to determine which fit is 
			better. For data that has prob_ct different concentration/
			probability pairs and niter different iterations of bootstrapped 
			probabilities:
		b_input: 2D np.array of size niters * 3. Should initially contain the
			guesses for b, and will ultimately contain the optimized values.
		prob_array: 2D np.array of size niters * prob_ct containing the 
			bootstrapped probability values.
		conc_list: 1D np.array of length prob_ct. Contains the concentrations
			that are matched to the probabilities. 
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b[0] and b[1]. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b[2].
		optim_method: An int from {0,1,2,3,4,5,6,7} to choose the solver. See
			the README for more details. 

		returns: 2D np.array of size niters * 3 containing the optimized 
			values for each bootstrap. When the ll2 curve is better, b[2] is 
			set to 1.
		'''
		if sigma_squared is None: sigma_squared = self.SS
		if beta_param is None: beta_param = self.BP
		if optim_method is None: optim_method = self.OPTIM_METHOD
		if lb is None: lb = self.LB
		if ub is None: ub = self.UB

		niters, prob_ct = prob_array.shape
		if self.use_C_lib: # pass off to the C function if it exists.
			funmin = np.zeros(niters)
			self.__ll23aAIC(prob_ct, niters, b_input, prob_array, conc_list,
						funmin, lb, ub, sigma_squared, 
						beta_param, optim_method)
		else:
			for i in range(niters):
				b_input[i] = self.min_llAIC(b_input[i], prob_array[i], 
									conc_list, sigma_squared, beta_param)
		return b_input

	def array_ls(self, nparam, b_input, prob_array, conc_list, 
										ls_method = None):
		'''
		Used to pass an array of probabilities to the C-based or python-based 
			least-squares solver. For data that has prob_ct different 
			concentration/probability pairs and niter different iterations of
			bootstrapped probabilities:
		nparam: Either 2 or 3, for whether the 2-paramter or 3-parameter dose-
			response curve should be used. 
		b_input: 2D np.array of size niters * 3 (yes, this must be this size
			even when using two parameters). Should initially contain the
			guesses for b, and will ultimately contain the optimized values.
		prob_array: 2D np.array of size niters * prob_ct containing the 
			bootstrapped probability values.
		conc_list: 1D np.array of length prob_ct. Contains the concentrations
			that are matched to the probabilities. 
		ls_method: An int from {0,1,2,3,4} to choose the solver. See the 
			README for more details. 

		returns: 2D np.array of size niters * 3 containing the optimized 
			values for each bootstrap.
		'''
		if ls_method is None: ls_method = self.LS_METHOD
		if self.use_C_lib: # pass off to the C function if it exists.
			niters, prob_ct = prob_array.shape
			self.__lsarray(nparam, prob_ct, niters, b_input, 
					conc_list, prob_array, ls_method)
		else:
			for i in range(niters):
				b_input[i] = FunctionFit.ll_ls(b_input[i], nparam, 
											prob_array[i], conc_list)
		return b_input

	@utils.surpress_warnings
	def min_ll2(self, b, probs, conc, sigma_squared = None, 
							optim_method = None, lb = None, ub = None):
		'''
		Used to pass nprob pairs of concentrations/probabilities to the 
			C-based or python-based ll2 optimizer.
		b: np.array of length 2. Should contain intial guess for
			optimal value.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b[0] and b[1]. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b[2].
		optim_method: An int from {0,1,2,3,4} to choose the solver. See the 
			README for more details. 

		returns: np.array of length 2 containing the optimized parameters.
		'''
		if sigma_squared is None: sigma_squared = self.SS
		if optim_method is None: optim_method = self.OPTIM_METHOD
		if lb is None: lb = self.LB
		if ub is None: ub = self.UB

		if self.use_C_lib:
			funmin = 0
			self.__ll2c(b, probs, conc, len(probs), funmin, lb[0:2], ub[0:2],
					sigma_squared, optim_method)
			return b
		else:
			bb = self.get_python_bounds(nparam = 2)
			res = minimize(FunctionFit.ll2p, b, 
					args = (probs, conc, sigma_squared), 
					method = 'SLSQP', jac = FunctionFit.ll2p_jac,
					bounds = bb)
			return res.x

	@utils.surpress_warnings
	def min_ll3(self, b, probs, conc, sigma_squared = None, 
					beta_param=None, optim_method = None, 
								lb = None, ub = None):
		'''
		Used to pass nprob pairs of concentrations/probabilities to the 
			C-based or python-based ll3 optimizer.
		b: np.array of length 3. Should contain intial guess for
			optimal value.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b[0] and b[1]. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b[2].
		optim_method: An int from {0,1,2,3,4} to choose the solver. See the 
			README for more details. 

		returns: np.array of length 3 containing the optimized parameters.
		'''

		if sigma_squared is None: sigma_squared = self.SS
		if beta_param is None: beta_param = self.BP
		if optim_method is None: optim_method = self.OPTIM_METHOD
		if lb is None: lb = self.LB
		if ub is None: ub = self.UB

		if self.use_C_lib:
			funmin = 0
			self.__ll3c(b, probs, conc, len(probs), funmin, lb, ub, sigma_squared,
					beta_param, optim_method)
			return b
		else:
			bb = self.get_python_bounds()
			res = minimize(FunctionFit.ll3p, b, 
					args = (probs, conc, sigma_squared, beta_param), 
					method = 'SLSQP', jac = FunctionFit.ll3p_jac,
					bounds = bb)
			return res.x

	def min_llAIC(self, b, probs, conc, sigma_squared = None, beta_param=None,
					optim_method = None, lb = None, ub = None):
		'''
		Used to pass nprob pairs of concentrations/probabilities to the 
			C-based or python-based ll2/ll3 optimizers and use the AIC to
			determine the optimal fit.
		b: np.array of length 3. Should contain intial guess for
			optimal value.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b[0] and b[1]. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b[2].
		optim_method: An int from {0,1,2,3,4} to choose the solver. See the 
			README for more details. 

		returns: np.array of length 3 containing the optimized parameters.
		'''
		if sigma_squared is None: sigma_squared = self.SS
		if beta_param is None: beta_param = self.BP
		if optim_method is None: optim_method = self.OPTIM_METHOD
		if lb is None: lb = self.LB
		if ub is None: ub = self.UB

		if self.use_C_lib:
			funmin = 0
			self.__ll23cAIC(b, probs, conc, len(probs), funmin, lb, ub,
					sigma_squared, beta_param, optim_method)
			return b
		else:
			bb = self.get_python_bounds()
			return self.ll23AIC_min(b, probs, conc, sigma_squared, 
				beta_param)

	@utils.surpress_warnings
	def min_ls(self, b, nparam, probs, conc, ls_method = None):
		'''
		Used to pass nprob pairs of concentrations/probabilities to the 
			C-based or python-based least-squares solver.
		b: np.array of length nparam. Should contain intial guess for
			optimal value.
		nparam: Either 2 or 3, for whether the 2-paramter or 3-parameter dose-
			response curve should be used. 
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		ls_method: An int from {0,1,2,3,4} to choose the solver. See the 
			README for more details. 

		returns:np.array of length nparam containing the optimized parameters.
		'''
		if ls_method is None: ls_method = self.LS_METHOD
		if self.use_C_lib: # pass off to the C function if it exists.
			nprob = len(probs)
			self.__lsfit(nparam, nprob, b, conc, probs, ls_method)
			return b
		else:
			return FunctionFit.ll_ls(b, nparam, probs, conc)

	@utils.surpress_warnings
	def min_cls(self, b, probs, conc, lb = None, ub = None):
		'''
		Used to pass nprob pairs of concentrations/probabilities to the 
			C-based or python-based constrained least-squares optimizer.
		b: np.array of length 3. Should contain intial guess for
			optimal value.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		optim_method: An int from {0,1,2,3,4} to choose the solver. See the 
			README for more details. 

		returns: np.array of length 2 containing the optimized parameters.
		'''
		if lb is None: lb = self.LB
		if ub is None: ub = self.UB

		if self.use_C_lib:
			funmin = 0
			self.__cls(b, probs, conc, len(probs), funmin, 
						lb, ub)
			return b
		else:
			bb = self.get_python_bounds()
			res = minimize(self.ls_SSR, b, args = (1.0-probs, conc), 
					method = 'SLSQP', bounds = bb)
			return res.x

	def array_cls(self, b_input, prob_array, conc_list, lb = None, ub = None):
		'''
		Used to pass an array of probabilities to the C-based or python-based
			ll3 solver. For data that has prob_ct different concentration/
			probability pairs and niter different iterations of bootstrapped 
			probabilities:
		b_input: 2D np.array of size niters * 3. Should initially contain the
			guesses for b, and will ultimately contain the optimized values.
		prob_array: 2D np.array of size niters * prob_ct containing the 
			bootstrapped probability values.
		conc_list: 1D np.array of length prob_ct. Contains the concentrations
			that are matched to the probabilities. 
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b[0] and b[1]. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b[2].
		optim_method: An int from {0,1,2,3,4,5,6,7} to choose the solver. See
			the README for more details. 

		returns: 2D np.array of size niters * 3 containing the optimized 
			values for each bootstrap. 
		'''
		if lb is None: lb = self.LB
		if ub is None: ub = self.UB

		niters, prob_ct = prob_array.shape
		if self.use_C_lib: # pass off to the C function if it exists.
			#funmin is for storing the function minimum, and can be returned
			#to provide the minimum value if needed. 
			funmin = np.zeros(niters)
			self.__cls_array(prob_ct, niters, b_input, prob_array, conc_list, 
						funmin, lb, ub)
		else:
			for i in range(niters):
				b_input[i] = self.min_cls(b_input[i], prob_array[i], conc_list)
		return b_input

	def ls_SSR(self, b, probs, conc):
		resid = self.least_squares_error(b, 3, probs, conc)
		return sum(resid*resid)

	def get_python_bounds(self, nparam = 3, LB = None, UB = None):
		'''
		Converts C-style bounds (vector of min and vector of max) to
		scipy bounds, which is a tuple of tuples of form
		((min,max), ..., (min,max))
		'''
		if LB is None: LB = self.LB
		if UB is None: UB = self.UB

		bb0 = (LB[0],UB[0])
		bb1 = (LB[1],UB[1])
		if nparam == 3:
			bb2 = (LB[2],UB[2])
			return (bb0, bb1, bb2)
		else:
			return (bb0, bb1)

	@staticmethod
	@utils.surpress_warnings
	def ll_ls(b, nparam, probs, conc):
		'''
		Fits the nparam-parameter dose-response curve (nparam = 2 or 3) using
			the method of least squares for the function
							       b2
						y = ------------------,
							1 + exp(b0 + b1*x)
			where b2=1.0 is fixed in the two-parameter, and free in the 3-
			parameter curves. 
			WARNING: For three parameters, b2 is not confined to the interval
			[0,1], and may therefore provide curves that are not statistically
			relevant. Use at your own risk. 
		b: np.array of length nparam. Should contain intial guess for
			optimal value.
		nparam: Either 2 or 3, for whether the 2-paramter or 3-parameter dose-
			response curve should be used. 
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.

		returns:np.array of length nparam containing the optimized parameters.
		'''
		if nparam == 2:
			return least_squares(FunctionFit.least_squares_error,
								b, args=(nparam, probs, conc)).x
		else:
			return least_squares(FunctionFit.least_squares_error,
								b, args=(nparam, probs, conc)).x

	@staticmethod
	@utils.surpress_warnings
	def ll2p(b, probs, conc, sigma_squared = SS):
		'''
		Log-likelihood function of the two parameter dose-response curve 
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein prior is b0, b1 ~ MVN(0, sigma*I2).
		
		b: np.array of length nparam. Should contain intial guess for
			optimal value.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b0 and b1. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b2.

		returns: the negative of the log-likelihood function (for optimization
		based on minimiation).
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
	def ll3p(b, probs, conc, sigma_squared = SS, beta_param=BP):
		'''
		Log-likelihood function of the three parameter dose-response curve 
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
		b2 ~ Beta(beta_param).

		b: np.array of length nparam. Should contain intial guess for
			optimal value.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b0 and b1. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b2.

		returns: the negative of the log-likelihood function (for optimization
		based on minimiation).
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
		'''
		Returns the AIC value of the log-likelihood function of both the 
		two- and three- parameter dose-response curves
								    b2
						y = ------------------
							1 + exp(b0 + b1*x)
		b: np.array of length nparam. Should contain optimized parameters for
			the dose-response curve.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b0 and b1. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b2.

		returns: float value of the AIC for the curve. 
		'''
		b0, b1, b2 = b

		#MVN Prior
		ll = -(b0*b0 + b1*b1)/(2*sigma_squared) 

		if (abs(b2 - 1) < 1e-6): #ll2
			p_sum = sum(probs)
			p_conc_sum = sum(probs*conc)
			ll = -(b0*b0 + b1*b1)/(2*sigma_squared) + b0*p_sum + b1*p_conc_sum - \
					sum(np.log(1 + np.exp(b0+b1*conc)))
			return 4 - 2*ll #AIC for two paramters
		else: #ll3
			xi = np.exp(b0+b1*conc)
			alpha = 1+xi
			ll += (beta_param[0]-1)*math.log(b2) + \
						(beta_param[1]-1)*math.log(1-b2) #Beta Prior
			ll += sum(probs*np.log(alpha-b2) - np.log(alpha) + \
						(1-probs)*math.log(b2))
			return 6 - 2*ll #AIC for three paramters

	@utils.surpress_warnings
	def ll23AIC_min(self, b, probs, conc, sigma_squared = SS, beta_param=BP):
		'''
		Cacluates the maximum of the log-likelihood function of both the 
		two- and three- parameter dose-response curves
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein priors are b0, b1 ~ MVN(0, sigma*I2) and b2 ~ Beta(beta_param)
		(three parameters) or when b2=1.0 (two parameters). The AIC is used
		to determine the optimal curve, the parameters of which are then 
		returned. 

		b: np.array of length nparam. Should contain intial guess for
			optimal value.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b0 and b1. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b2.

		returns:np.array of length nparam containing the optimized parameters.
		'''
		#fit ll2
		bb = self.get_python_bounds(nparam = 2)
		res2 = minimize(FunctionFit.ll2p, b[0:2], 
					args = (probs, conc, sigma_squared), 
					method = 'SLSQP', jac = FunctionFit.ll2p_jac, bounds = bb)
		b2 = np.array([*res2.x, 1.0])
		AIC2 = FunctionFit.ll_for_AIC(b2, probs, conc, sigma_squared)
		#fit ll3
		bb = self.get_python_bounds(nparam = 3)
		res3 = minimize(FunctionFit.ll3p, b, 
					args = (probs, conc, sigma_squared, beta_param), 
					method = 'SLSQP', jac = FunctionFit.ll3p_jac, bounds = bb)
		AIC3 = FunctionFit.ll_for_AIC(res3.x, probs, conc, 
									sigma_squared, beta_param)

		if (AIC2 < AIC3):
			return b2
		else:
			return res3.x

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
		
		b: np.array of length nparam containing the paramter values at which
			to calculate the Jacobian.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b0 and b1. 

		returns: np.array of length nparam containing the optimized parameters.
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
	def ll3p_jac(b, probs, conc, sigma_squared = SS, beta_param=BP):
		'''
		Jacobian of the log-likelihood function of the three 
		parameter dose-response curve 
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
		b2 ~ Beta(beta_param).

		b: np.array of length nparam containing the paramter values at which
			to calculate the Jacobian.
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b0 and b1. 
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b2.

		returns: np.array of length nparam containing the optimized parameters.
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
	def least_squares_error(b, nparam, probs, conc):
		'''
		Dose-response curve for least-squares fitting. If len(b)=2, then 
		b2 = 1; else if len(b)=3, then b2 = b[2]
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		b: np.array of length nparam containing the parameter values at which
			to calculate the error.
		nparam: Either 2 or 3, for whether the 2-paramter or 3-parameter dose-
			response curve should be used. 
		probs: np.array of length nprob containing the probability values.
		conc: np.array of length nprob containing the concentration values.

		returns: float containing the summed squared residuals.
		'''
		if nparam == 2:
			return np.array((1-1/(1 + np.exp(b[0] + b[1]*conc)) - probs))
		else:
			if b[2] > 1:
				return 1e10
			else:
				return np.array((1-b[2]/(1.0 + np.exp(b[0] +\
												b[1]*conc)) - probs))

	@staticmethod
	@utils.surpress_warnings
	def loglogit2(b, conc): 
		'''
		Calculates the values of the two parameter dose-response curve 
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		given b = [b0, b1] and x.
		b: np.array of length 2 containing the parameter values at which
			to evaluate the function.
		conc: 1D np.array containing the concentration values at which to 
			evaluate the function.

		returns: 1D np.array of length len(conc) containing the function 
			values at each concentration 
		'''
		return 1./(1.+ np.exp(b[0] + conc * b[1]))

	@staticmethod
	@utils.surpress_warnings
	def loglogit3(b, conc): 
		'''
		Calculates the values of the three parameter dose-response curve 
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		given b = [b0, b1, b2] and x.

		b: np.array of length 3 containing the parameter values at which
			to evaluate the function.
		conc: 1D np.array containing the concentration values at which to 
			evaluate the function.

		returns: 1D np.array of length len(conc) containing the function 
			values at each concentration 
		'''
		return b[2]/(1.+ np.exp(b[0] + conc * b[1]))

	@staticmethod
	def estimate_initial_b(conc, probs, params = 3, rev = False):
		'''
		Produce an initial estimate of the starting vector for curve 
		optimization. The slope defualts to 1, and the data is used to 
		generate estimates of the baseline mortality and the LC50. 
		
		conc: 1D np.array containing the concentration values at which to 
			evaluate the function.
		probs: np.array of length len(conc) containing the probability values.
		params: Either 2 or 3; number of parameters in curve. 
		returns: np.array of len(params) for initial guess of values. 
		'''		
		#no great way to estimate slope implemented yet without a curve fit.
		default_slope= 1
		#estimate background mortality:
		n_vals = round(0.2 * len(conc))
		#sort the lists
		zipped =  sorted(zip(conc, probs))
		tuples = zip(*zipped)
		conc, probs = [list(val) for val in  tuples]
		conc = np.array(conc)
		probs = np.array(probs)

		#set the lc50 y-value with(out) background mortality in consideration.
		med = 0.5
		if params == 3:
			background_mort = sum( probs[:n_vals])/n_vals #ave
			if rev: background_mort = 1-background_mort
			med = background_mort/2.
			#prevent nasty numerical errors, e.g. undefined function values
			#from beta distributions when background_mort >=1 or <=0.
			if abs(background_mort - 1.0)<1e-6: background_mort = 0.99
			if background_mort < 0.5 : background_mort = 0.75

		#estimate the b0 parameter
		high_idx = np.where(probs > med)[0]

		if len(high_idx) == 0: 
			est_intercept = max(conc)-1
		else:
			est_pt = round(len(high_idx)*.6)-1
			est_intercept= - conc[high_idx[est_pt]]/default_slope

		if params == 3:
			return np.array([est_intercept, default_slope, background_mort])
		else:
			return np.array([est_intercept, default_slope])


	def switch_fitter(self, switch, conc, probs, 
				sigma_squared = None, beta_param = None,
				rev = True, perturbation = False):
		'''
		Drives the selection of the apprpriate function to fit based on the 
		value specified in the switch. 
		switch: selected from "auto", "ls2", "ls3", "ll2", "ll3", "AIC" in 
			anlysis_config.txt
		conc: The np.array of concentrations; shape = (nconc,)
		props: The np.array of survival propbabilities. shape = (nconc,)
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b0 and b1. To be passed to bayesian fittings
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b2. To be passed to bayesian fittings
		Returns: np.array of fitted parameters. shape = (3,)
		'''
		if sigma_squared is None: sigma_squared = self.SS
		if beta_param is None: beta_param = self.BP
		b2 = self.estimate_initial_b(conc, probs, params = 2, rev = rev)
		b3 = self.estimate_initial_b(conc, probs, params = 3, rev = rev)
		if perturbation: 
			b2 += np.random.normal(0, 0.5, 2)
			b3 += np.array([*np.random.normal(0, 0.5, 2) , 0.0])

		switch = switch.lower()

		if switch == 'auto':
			if b3[2] > 0.10: switch = "ll3"
			else: switch = "ll2"

		if switch in ["2", "ll2", 2]:
			b = self.min_ll2(b2, probs, conc, 
				sigma_squared = sigma_squared)
		elif switch in ["3", "ll3", 3]:
			b = self.min_ll3(b3, probs, conc,
				sigma_squared = sigma_squared, beta_param = beta_param)
		elif switch in ["ls3"]:
			b = self.min_ls(b3, 3, probs, conc)
		elif switch in ["ls2"]:
			b = self.min_ls(b2, 2, probs, conc)
		elif switch in ["cls", "cls3"]:
			b = self.min_cls(b3, 1.0 - probs, conc)
		elif switch in ["best", "aic"]:
			b = self.min_llAIC(b3, probs, conc,
				sigma_squared = sigma_squared, beta_param = beta_param)
		else: #set default behavior to LL3
			b = self.min_ll3(b3, probs, conc,
				sigma_squared = sigma_squared, beta_param = beta_param)

		if len(b) == 2:
			return np.array([b[0], b[1], 1.0])
		else:
			return b

	def switch_array_fitter(self, switch, conc, prob_array, init_prob,
				sigma_squared = None, beta_param = None,
				perturbation = False):
		'''
		Drives the selection of the apprpriate function to fit based on the 
		value specified in the switch. 
		switch: selected from "auto", "ls2", "ls3", "ll2", "ll3", "AIC" in 
			anlysis_config.txt
		conc: The np.array of concentrations; shape=(nconc, )
		prob_array: The np.array of survival propbabilities; shape=(n_iter,nconc)
		init_prob: The np.array of true probs for estimating the initial value 
			of b to start optimization; shape=(nconc,)
		sigma_squared: float value for sigma**2 for the MVN prior for the
			values of b0 and b1. To be passed to bayesian fittings
		beta_param: 1D np.array of length 2 contain the two beta parameters
			for the prior on b2. To be passed to bayesian fittings
		Returns: (n_iter,3) np.array of fitted parameters. 
		'''
		if sigma_squared is None: sigma_squared = self.SS
		if beta_param is None: beta_param = self.BP

		#generate initial estimates
		niters, nprobs = prob_array.shape
		b2 = self.estimate_initial_b(conc, init_prob, params = 2, rev = True)
		b2_array = np.repeat([b2], niters, axis=0)
		b3 = self.estimate_initial_b(conc, init_prob, params = 3, rev = True)
		background_mort = b3[2]
		b3_array = np.repeat([b3], niters, axis=0)

		if perturbation:
			b2_array += np.random.normal(0, 0.5, b2_array.shape)
			b3_array += np.c_[np.random.normal(0, 0.5, b2_array.shape), 
								np.zeros(niters)]

		switch = switch.lower() #verify switch is lower-case for matching

		if switch == 'auto':
			if background_mort > 0.10: switch = "ll3"
			else: switch = "ll2"

		if switch in ["2", "ll2", 2]:
			b_out = self.array_ll2(b2_array, prob_array, conc, 
				sigma_squared = sigma_squared)
		elif switch in ["3", "ll3", 3]:
			b_out = self.array_ll3(b3_array, prob_array, conc,
				sigma_squared = sigma_squared, beta_param = beta_param)
		elif switch in ["ls3"]:
			b_out = self.array_ls(3, b3_array, prob_array, conc)
		elif switch in ["ls2"]:
			#NB: note that whilst two parameters are used, the array still
			#is required to be three param wide. 
			b_out = self.array_ls(2, b3_array, prob_array, conc)
		elif switch in ["cls", "cls3"]:
			b_out = self.array_cls(b3_array, 1.0 - prob_array, conc)
		elif switch in ["best", "aic"]:
			b_out = self.array_ll23AIC(b3_array, prob_array, conc,
				sigma_squared = sigma_squared, beta_param = beta_param)
		else:
			b_out = self.array_ll3(b3_array, prob_array, conc,
				sigma_squared = sigma_squared, beta_param = beta_param)

		biter, blen = b_out.shape
		#need to make an niters * 3 array
		if blen == 2:
			return np.c_[b_out,np.ones(biter)]
		else:
			return b_out