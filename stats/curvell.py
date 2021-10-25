# stats/curvell.py
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
from scipy.optimize import minimize, least_squares
from scipy.interpolate import CubicSpline
import numpy as np
from stats.corr_beta import corr_beta_vars
import matplotlib.pyplot as plt
import warnings
from statistics import median, mean
from sklearn.neighbors import KernelDensity
from stats.merlin_grapher import MerlinGrapher
from time import time
import stats.utils as utils
from itertools import compress
from stats.functionfit import FunctionFit

def default_params_dict():
	param_dict = {"BOOTSTRAP_ITERS": 1000,  
						"CURVE_TYPE": 'best', 
						"FIT_METHOD": 'BFGS', 
						"BETA_PRIOR": 0.1, 
						"RHO": 0.1,
						"N_POINTS": 151
				 }
	return param_dict

class CI_finder:	
	'''
	Performs the curve-fitting and associated math for dose-response analysis
	for the bioassay.
	'''
	default_options = default_params_dict()

	def __init__(self, 
				unique_plate_ids, 
				live_count, 
				dead_count, 
				conc, 
				include_now, 
				options = None, 
				**kwargs):
		if unique_plate_ids is None or live_count is None or \
				dead_count is None or conc is None:
			raise ValueError("CI_finder requries plate_ids, live_count, " +\
				"dead_count, and conc values.")
		elif not (len(unique_plate_ids)== len(live_count) and 
			len(unique_plate_ids)==len(dead_count) and 
			len(unique_plate_ids)==len(conc) and 
			len(unique_plate_ids)==len(include_now)):
			raise ValueError("CI_finder requries plate_ids, live_count, " +\
				"dead_count, and conc to all be the same length.")
		self.unique_plate_ids =  list(compress(unique_plate_ids, include_now))
		self.live_count = np.array(list(compress(live_count, include_now)))
		self.dead_count = np.array(list(compress(dead_count, include_now)))
		self.conc = np.array(list(compress(conc, include_now)))
		
		#for the purpose of determining if there are new UIDs included
		self.uid_list = list(set(self.unique_plate_ids))

		self.n_trials = kwargs['n_trials']
		self.options = self.default_options
		if options:
			for k, v in options.items(): self.options[k] = v
		#values created by calculation
		self.params = None
		self.points = None
		self.plot_quant = None

	#DECORATOR
	def update_options(func, *args, **kwargs): 
		'''
		Decorator to update the options for this class if the options
		keyword is included. 
		'''
		def inner(*args, **kwargs):
			if "options" in kwargs:
				for k, v in options.items(): self.options[k] = v
			func(*args, **kwargs)
		return inner

	@staticmethod
	@utils.surpress_warnings
	def ll3(b, probs, conc, sigma_squared = 1e6, beta_param=[1.5,1.01]):
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
	def least_squares_fit(b, probs, conc):
		'''
		Dose-response curve for least-squares fitting. If len(b)=2, then 
		b2 = 1; else if len(b)=3, then b2 = b[2]
							    b2
					y = ------------------
						1 + exp(b0 + b1*x)
		'''
		if len(b) == 2:
			return np.array((1-1/(1 + np.exp(b[0] + b[1]*conc)) - probs)**2)
		else:
			if b[2] > 1:
				return 1e10
			else:
				return np.array((1-b[2]/(1 + np.exp(b[0] + b[1]*conc)) - probs)**2)


	@staticmethod
	@utils.surpress_warnings
	def ll3_jac(b, probs, conc,sigma_squared = 1e6, weibull_param=[2,1]):
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
	def ll2(b, probs, conc, sigma_squared = 1e6):
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
	def ll2_jac(b, probs, conc, sigma_squared = 1e6):
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

	@staticmethod
	@utils.surpress_warnings
	def loglogit2(b, conc): 
		'''
		Calculates the values of the two parameter dose-response curve 
							    1
					y = ------------------
						1 + exp(b0 + b1*x)
		given b = [b0, b1] and x.
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
		'''
		return b[2]/(1.+ np.exp(b[0] + conc * b[1]))

	def ll3_find_LC(self, quantiles):
		'''
		Given the list self.params of parameters from bootstrapping and the 
		list of quantiles, this method returns the concentrations at which each
		of the quantiles is met. 
		'''
		quant2 = np.tile(np.array(quantiles), (len(self.params),1))
		params = np.reshape(np.repeat(np.array(self.params[:,0:2]), 
					len(quantiles)), 
					(self.params[:,0:2].shape[0], 
							self.params[:,0:2].shape[1], 
							len(quantiles)))
		return (np.log(1./quant2 - 1.)-params[:,0])/params[:,1]

	@staticmethod
	def estimate_initial_b(conc, probs, params = 3, rev = False):
		'''
		Produce an initial estimate of the starting vector for curve 
		optimization. The slope defualts to 1, and the data is used to 
		generate estimates of the baseline mortality and the LC50. 
		'''
		#no good way to estimate slope yet without a curve fit.
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

		#estimate the b0 parameter
		high_idx = np.where(probs > med)[0]

		if len(high_idx) == 0: 
			est_intercept = max(conc)
		else:
			est_intercept= -conc[high_idx[-1]]/default_slope

		if params == 3:
			return np.array([est_intercept, default_slope, background_mort])
		else:
			return np.array([est_intercept, default_slope])

	def fit_curve(self, probs):
		'''
		Curve-fitting driver. 
		'''
		#TODO: better methods of b estimation
		#TODO: handle cases where EC50 is not going to be in the data :( 
		
		@utils.surpress_warnings
		def fit_ll3(b, probs, functionFitter):
			'''
			Inner function that attempts to fit the three-paremeter dose-
			response curve with the option-specified method for minimization.
			If this method fails, the more computationally-expensive but 
			better-behaved Nelder-Mead method is used. 
			'''

			# res = minimize(self.ll3, b, args = (probs, self.conc), 
			# 		method = self.options["FIT_METHOD"], jac = self.ll3_jac)
			# if not res.success:
			# 	res = minimize(self.ll3, b, args = (probs, self.conc), 
			# 		method = 'Nelder-Mead')
			# return res
			return functionFitter.min_ll3(b, probs, self.conc )
		
		functionFitter = FunctionFit()
		b2 = self.estimate_initial_b(self.conc, probs, params = 2, rev = True)
		b3 = self.estimate_initial_b(self.conc, probs, params = 3, rev = True)

		if self.options["CURVE_TYPE"].lower() == 'auto':
			if background_mort > 0.10: self.options["CURVE_TYPE"] = "ll3"
			else: self.options["CURVE_TYPE"] = "ll2"

		if self.options["CURVE_TYPE"].lower() in ["2", "ll2", 2]:
			return minimize(self.ll2, b2, args = (probs, self.conc), 
					method = self.options["FIT_METHOD"], jac = self.ll2_jac).x
		elif self.options["CURVE_TYPE"].lower() in ["3", "ll3", 3]:
			return fit_ll3(b3, probs,functionFitter)
		elif self.options["CURVE_TYPE"].lower() in ["ls3"]:
			return least_squares(self.least_squares_fit, b3, args=(probs, self.conc)).x
		elif self.options["CURVE_TYPE"].lower() in ["ls2"]:
			return least_squares(self.least_squares_fit, b2, args=(probs, self.conc)).x
		elif self.options["CURVE_TYPE"].lower() in ["best", "aic"]:
			res2 = minimize(self.ll2, b2, args = (probs, self.conc), 
					method = self.options["FIT_METHOD"], jac = self.ll2_jac)
			b3p = np.array([res2.x[0], res2.x[1], background_mort])
			res3 = fit_ll3(b3p, probs, functionFitter)
			AIC2 = 4 - 2*res2.fun
			AIC3 = 6 - 2*res3.fun
			return res2.x if AIC2 <= AIC3 else res3.x

	def bootstrap_CIs(self):
		'''
		Driver for the bootstrapping process. Because the curvefitting is 
		computationally costly but is embrassingly parallel, this is best 
		done in parallel. 
		'''
		import multiprocessing
		from joblib import Parallel, delayed

		if self.params is not None: return
		#get correlated beta variables
		beta_probs = corr_beta_vars(self.unique_plate_ids, 
									self.live_count, 
									self.dead_count, 
									size = self.options["BOOTSTRAP_ITERS"], 
									scale = self.options["BETA_PRIOR"], 
									null_scale = self.options["BETA_PRIOR_0"],
									rho = self.options["RHO"])
		#calculate point data points while we have the beta parameters. 
		self.get_point_error_bars(beta_probs)
		#possibility of having 2 or 3 parameters
		self.params = np.zeros((self.options["BOOTSTRAP_ITERS"], 3))
		cpu_count = multiprocessing.cpu_count() 
		if self.options["BOOTSTRAP_ITERS"] > cpu_count : 
			self.options["BOOTSTRAP_ITERS"] = cpu_count
		dict_list = Parallel(n_jobs = self.options["BOOTSTRAP_ITERS"])(delayed(
					self.fit_curve)(row) for row in beta_probs)
		# dict_list = Parallel(n_jobs = 1)(delayed(
		# 			self.fit_curve)(row) for row in beta_probs)
		for iter_count, res in enumerate(dict_list):
			if len(res) == 3:
				self.params[iter_count] = res
			else:
				self.params[iter_count,0:2] = res
				self.params[iter_count,2] = 1.
		
	def get_point_error_bars(self, beta_probs):
		'''
		Calculates the error bars for each data point based on the beta
		variables from the bootstrap samples. 
		'''
		lb = (1 - self.options["ERROR_BAR_CI"])/2.
		ub = 1-lb
		errors = np.quantile(beta_probs, [ub, lb], 
				interpolation='linear',
				axis = 0)
		probs = self.dead_count / (self.dead_count+self.live_count)
		probs = np.tile(np.array(probs), (2,1))
		self.error_bars = np.abs(errors-probs)

	def calculate_curves(self):
		'''
		Calculates the value of the dose response curves at each value in
		self.x for each set of parameters. 
		'''
		if self.points is not None: return
		self.min_conc, self.max_conc = min(self.conc)-1, max(self.conc)+1
		lb, ub = math.floor(self.min_conc), math.ceil(self.max_conc)
		self.x = np.linspace(lb-1, ub+1, self.options["N_POINTS"])
		self.points = np.zeros((len(self.params), self.options["N_POINTS"]))
		for iter_count, row in enumerate(self.params): 
			self.points[iter_count] = self.loglogit3(row, self.x)
		self.r2 = self.find_r2(fit_x = self.x, 
			fit_y = np.median(self.points, axis = 0), 
			conc = self.conc, 
			probs = self.live_count / (self.dead_count+self.live_count))
		# print(self.r2)

	@staticmethod
	def find_r2(fit_x, fit_y, conc, probs):
		'''
		Finds the r**2 value of the curve using the formula 
		r**2 = 1 - SS_{res}/SS_{tot}, where SS_{res} = sum(y_i - f_i)^2 for 
		y_i being the response at concentration i and the f_i the value of the
		dose-response at i, and SS_{tot} = sum(y_i - y_bar)^2, with y_bar 
		being the average response across all concentrations. 
		'''
		spline = CubicSpline(fit_x, fit_y)
		spline_vals = spline(conc)
		start = 0
		end = 10
		SSR = spline_vals - probs
		sum_of_square_resids = sum(SSR * SSR)
		SSF = probs - mean(probs)
		sum_of_square_nofit = sum(SSF*SSF)
		return 1- sum_of_square_resids/sum_of_square_nofit

	@update_options
	def get_plot_CIs(self, quantiles = [.025, 0.5, 0.975], options=None):
		'''
		Driver for calculating the confidence intervals for the plot.
		'''
		if self.plot_quant is None: 
			self.bootstrap_CIs()
			self.calculate_curves()
		self.plot_quant = utils.calc_ET_CI(self.points, 
				CI_level = self.options["LC_CI"], 
				resample = False)

	def get_CIs(self, LC_VALUES = None, LC_CI = 0.95, log = False):
		'''
		Driver for calculating dose-response intervals. 
		'''

		if self.params is None: self.bootstrap_CIs()
		if LC_VALUES is None: LC_VALUES = 1-self.options["LC_VALUES"]
		if LC_CI is None: LC_CI = self.options["LC_CI"]
		EC_vals = self.ll3_find_LC(LC_VALUES)
		EC_CI = self.get_EC_CIs(EC_vals, LC_CI)
		# print("LC_CI:", LC_CI)
		# print("LC_VALUES:", LC_VALUES)
		# print("LC_VALUES:", np.power(2.,EC_CI))
		return EC_CI if log else np.power(2.,EC_CI)

	def get_EC_CIs(self, EC_vals, CI_val):
		'''
		Finds the confidence intervals for doses. 
		'''
		# print(self.options["CI_METHOD"])
		func = utils.calc_ET_CI if self.options["CI_METHOD"].lower() in \
				utils.ET_VARS else utils.calc_HPDI_CI
		EC_val_summary = func(EC_vals, CI_level = CI_val)
		return EC_val_summary.T if EC_val_summary.ndim>1 else EC_val_summary

	#Returns the LC50 credible interval 
	def get_LC50_CI(self, CI_val=0.95, log = False): 
		return  self.get_CIs(LC_VALUES = np.array([0.5]), 
			LC_CI = 0.95, log=log).squeeze()

	#Returns a slope credible interval via self.get_param_CI.
	def get_slope_CI(self, CI_val=0.95): 
		return self.get_param_CI(1, CI_val).reshape((-1))

	#Returns the baseline moretality estimate via self.get_param_CI.
	def get_baseline_mort_CI(self, CI_val=0.95): 
		return self.get_param_CI(2, CI_val).reshape((-1))

	def get_param_CI(self, parameter, CI_val):
		'''
		Produces a confidence interval for a parameter based on the 
		bootstrapped fits. 
		'''
		if self.params is None: self.bootstrap_CIs()
		func = utils.calc_ET_CI if self.options["CI_METHOD"].lower() in \
				utils.ET_VARS else utils.calc_HPDI_CI
		return func(self.params[:,parameter], CI_level = CI_val)

	def LC_kernel(self, LC_val = 0.5):
		'''
		Generates a KernelDensity object from an LC value.
		'''
		vals = self.ll3_find_LC(quantiles = LC_val)
		return KernelDensity(kernel ='gaussian', bandwidth = 0.5).fit(vals)

	def plot_CIs(self):
		'''
		Produces a MerlinGrapher object set up with options and data for 
		plotted points, error bars, best fit line, and line CI. 
		'''
		lb = (1 - self.options["CURVE_CI"])/2
		ub = 1 - lb
		self.get_plot_CIs(quantiles = [lb, 0.5, ub])
		merlin_plot = MerlinGrapher(x = self.x, 
									lb = self.plot_quant[0],
									ub = self.plot_quant[2],
									line = self.plot_quant[1],
									conc = self.conc,
									probs = self.live_count/(self.dead_count +\
											self.live_count), 
									error_bars = self.error_bars,
									n_trials = self.n_trials,
									options = self.options)
		return merlin_plot

	def reset_curves(self):
		'''
		This function clears out the params, points, and CI lines to allow one
		do new bootstrapping.
		'''
		self.params = None
		self.points = None
		self.plot_quant = None
