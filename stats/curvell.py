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
						"N_POINTS": 151,
						'LS_METHOD': 4,
						'OPTIM_METHOD': 5,
						'LL_SIGMA': 1000,
						'LL_BETA1': 1.5,
						'LL_BETA2': 1.001,
						'USE_CLIB': True
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
		# print(options.items())
		if options is not None:
			for k, v in options.items(): 
				self.options[k] = v
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

	def ll3_find_LC(self, quantiles):
		'''
		Given the list self.params of parameters from bootstrapping and the 
			list of quantiles, this method returns the concentrations at which 
			each of the quantiles is met. 
		quantiles: A list-like of quanties for calculating concentrations.
		returns: 1D np.array of concentrations at which each concentraion is
			met.
		'''
		quant2 = np.tile(np.array(quantiles), (len(self.params),1))
		params = np.reshape(np.repeat(np.array(self.params[:,0:2]), 
					len(quantiles)), 
					(self.params[:,0:2].shape[0], 
							self.params[:,0:2].shape[1], 
							len(quantiles)))
		return (np.log(1./quant2 - 1.)-params[:,0])/params[:,1]

	def fit_curve(self, probs):
		'''
		Wrapper for fitting a single curve to a single set of concentrations/
			probabilities. Curve type is specified via analysis_config.txt. 
			Fitting may occur in either a C library or in Python. 
			!!WARNING!!: in Windows OS, there is the potential to have serious
			function referencing issues, as function references may be capped
			at 1024 or 2048, preventing the number of bootstraps from being 
			higher than this!!
		probs: 1D np.array of probabilities to fit. 
		returns: np.array of length 3 containing the paramters of the fitted
			dose response curve.
		'''
		ff = FunctionFit(**self.options)
		switch = self.options["CURVE_TYPE"].lower()
		return ff.switch_fitter(switch, self.conc, probs)

	def array_curve(self, beta_probs):
		'''
		Wrapper for fitting multiple curves to a multiple sets of 
			concentration/probability pairs. Curve type is specified via 
			analysis_config.txt. Fitting may occur in either a C library or in
			Python. Performing this fitting in arrays is crucial to avoid 
			referencing issues on Windows systems where the number of 
			references to C functions is limited to 1024 or 2048, which may be
			well below the number of bootstraps requested. 
		beta_probs: 2D np.array of shape niters*nprobs of probabilities to 
			fit, with niters being the number of different sets of 
			probabilities, and nprobs being the number of conc/prob pairs.
		returns: 2D np.array of shape (niters,3) containing the paramters of  
			the fitted dose response curves.
		'''

		ff = FunctionFit(**self.options)
		#init_prob is included for calculating the initial vector of parameters
		#for curve fitting.
		init_prob = self.dead_count / (self.dead_count+self.live_count)
		switch = self.options["CURVE_TYPE"].lower()
		return ff.switch_array_fitter(switch, self.conc, beta_probs, init_prob)

	def bootstrap_CIs(self):
		'''
		Driver for the bootstrapping process. Because the curvefitting is 
			computationally costly but is embrassingly parallel, this is best 
			done in parallel. 
		The np.array of parameters for the fitted curve is of shape 
			 (n_bootstraps, 3), and is saved in attribute self.params
		'''
		import multiprocessing
		from joblib import Parallel, delayed

		niters = self.options["BOOTSTRAP_ITERS"]

		if self.params is not None: return #already calcuated. 

		#First create an array of partially-correlated beta variables, as 
		#defined in the statistics summary .pdf document. 
		beta_probs = corr_beta_vars(self.unique_plate_ids, 
									self.live_count, 
									self.dead_count, 
									size = niters, 
									scale = self.options["BETA_PRIOR"], 
									null_scale = self.options["BETA_PRIOR_0"],
									rho = self.options["RHO"])
		#calculate point data points whilst we have the beta parameters. 
		self.get_point_error_bars(beta_probs)
		cc = multiprocessing.cpu_count() 
		#if user wants more cpus than exist, set to the number of actual cpus
		cu = cc if self.options["NCPU"] > cc else self.options["NCPU"]
		if cu <=0: cu += cc + 1 #in case cpus is set to a negative number 
				
		#create a FunctionFit object, and then set the default paramters for
		#bayesian curve fitting for use as needed. 
		ff = FunctionFit(**self.options)
		ff.SS = self.options["LL_SIGMA"] * self.options["LL_SIGMA"]
		ff.BP=np.array([self.options["LL_BETA1"],self.options["LL_BETA2"]])

		#mainly for debugging purposes; setting this
		use_array_method = True  

		if use_array_method:
			#first make a list of arrays for parallel computing
			# make divisions based on number of cpus to use
			groups = cu*2
			if groups > niters: #there are more iters
				#than there are bootstrap iters... 
				groups = 0;
				#so set the number of cpus to 1. 
				cu = 1
			'''
			To distribute all bootstrapped sets of probabilities to the 
			different processors, we make a list of different 2D np.arrays of
			the probabilities, and use this in Parallel. 
			'''
			breaks = np.round(np.linspace(0,niters,groups+1))
			beta_prob_list = []
			for i in range(groups):
				beta_prob_list.append(beta_probs[int(breaks[i]):\
														int(breaks[i+1]),:])
			dict_list = Parallel(n_jobs = cu)(delayed(self.array_curve)\
								(beta_probs) for beta_probs in beta_prob_list)
			#unpack the dictionary, which involves converting a 
			#list of list of answers to np.array of answers
			self.params =  np.array([item for answer in dict_list for \
															item in answer])

		else:
			#possibility of having 2 or 3 parameters
			self.params = np.zeros((niters, 3))
			dict_list = Parallel(n_jobs = cu)(delayed(
						self.fit_curve)(row) for row in beta_probs)
			#unpack the dictionary
			for iter_count, res in enumerate(dict_list):
				if len(res) == 3:
					self.params[iter_count] = res
				else:
					self.params[iter_count,0:2] = res
					self.params[iter_count,2] = 1.

	def get_point_error_bars(self, beta_probs):
		'''
		Calculates the error bars for each data point based on the beta
		variables from the bootstrap samples. The width of the error bars is
		set by ERROR_BAR_CI in analysis_config.txt, and the error bars are
		saved in the attribute self.error_bars.
		beta_probs: a 2D array of shape (bootstrap_iters, n_probs)
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
			self.points[iter_count] = FunctionFit.loglogit3(row, self.x)
		self.r2 = self.find_r2(fit_x = self.x, 
			fit_y = np.median(self.points, axis = 0), 
			conc = self.conc, 
			probs = self.live_count / (self.dead_count+self.live_count))

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
		Generates a KernelDensity object from an LC value. Estimates the
			bandwidth for the gaussian kernel density using 
			utils.est_gaussian_kernel_bandwidth
		LC_val: The LC value for which the kernel should be created. 
		returns: a KernelDensity object  
		'''
		vals = self.ll3_find_LC(quantiles = [LC_val])
		bw = utils.est_gaussian_kernel_bandwidth(vals)
		return KernelDensity(kernel ='gaussian', bandwidth = bw).fit(vals)

	def plot_CIs(self):
		'''
		Produces a MerlinGrapher object set up with options and data for 
		plotted points, error bars, best fit line, and line CI. 
		Returns a MerlinGrapher object intialized with all the data needed
		for generating the graphs that are used in the LaTeX output. 
		'''
		lb = (1 - self.options["CURVE_CI"])/2
		ub = 1 - lb
		self.get_plot_CIs(quantiles = [lb, 0.5, ub])
		merlin_plot = MerlinGrapher(x = self.x, 
									lb = self.plot_quant[0],
									ub = self.plot_quant[2],
									line = self.plot_quant[1],
									conc = self.conc,
									probs = self.live_count/(self.dead_count+\
											self.live_count), 
									error_bars = self.error_bars,
									n_trials = self.n_trials,
									options = self.options)
		return merlin_plot

	def reset_curves(self):
		'''
		This function clears out the params, points, and CI lines to allow one
		do new bootstrapping. This is usefule when data is added or removed 
		and curves need to be recalculated. 
		'''
		self.params = None
		self.points = None
		self.plot_quant = None
