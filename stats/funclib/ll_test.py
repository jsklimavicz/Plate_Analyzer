import numpy as np
import math
from ctypes import *
import time 
from scipy.optimize import minimize, least_squares


def ll3(b, probs, conc, sigma_squared = 1e6, beta_param=[1.5,1.01]):
	'''
	Log-likelihood function of the three parameter dose-response curve 
						    b2
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
	b2 ~ Weibull(weibull_params).

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

	ll2 = (b0*b0 + b1*b1)/(2*sigma_squared)
	# print(f"Before start: ll = {ll2}")

	# for i in range(len(conc)):
	# 	ll2 -= probs[i]*b0 + b1*probs[i]*conc[i] - math.log(1+math.exp(b0+b1*conc[i]))
	# 	print(f"After {i}: ll = {ll2}")

	return(-ll)

def ll3_jac(b, probs, conc,sigma_squared = 1e6, beta_param=[1.5,1.01]):
	'''
	Jacobian of the log-likelihood function of the three 
	parameter dose-response curve 
						    b2
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
	b2 ~ Weibull(weibull_params).
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

if __name__ == "__main__":

	so_file = "./cloglik.so"
	cloglik = CDLL(so_file)

	conc = np.array([8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125]*4)
	probs = np.array([0, 0, 0, .33, 0.6875, 0.846153846153846, 1, 0.857142857142857, 0.923076923076923, 0.914545454545455, 0, 0, 0, 0, 0, .2, .43, 0.583333333333333, 0.85, 0.916666666666667]*4)
	probs2 = np.sqrt(probs)
	probs3 = probs *0.8 + 0.01
	probs_full = np.array([probs,probs2,probs3])
	n_reps = 1
	probs_full = np.repeat(probs_full,n_reps, axis=0)
	# print(probs_full)
	b3 = [.89,-1.3,.917]
	b3_full = np.array([b3]*3)
	b3_full = np.repeat(b3_full,n_reps, axis=0)
	b3_empty = np.zeros_like(b3_full)

	b2 = [.89,-1.3]
	b2_full = np.array([b2]*3)
	b2_full = np.repeat(b2_full,n_reps, axis=0)
	b2_empty = np.zeros_like(b2_full)


	# n=10
	# t0 = time.time()
	# for i in range(n):
	# 	for j in range(len(probs_full)):
	# 		res = minimize(ll3, b_full[j], args = (probs_full[j], conc), method = 'BFGS', jac = ll3_jac)
	# 		b_empty[j] = res.x
	# # print(b_full)
	# t1 = time.time()
	# print(f"LL3 Numpy: {(t1-t0)/n*1e6:.3f} us")
	# # print(f"LL3 Numpy: {(t1-t0)/n*1e6:.3f} us; Grad: {ll3_jac(np.array(res.x), np.array(probs), np.array(conc))}")



	# ll3c = cloglik.ll3_min
	# b_type = c_double * len(b)
	# data_type = c_double * len(probs)
	# beta_type = c_double * 2
	# funmin = c_double(0)
	# ll3c.argtypes = (np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
	# 				np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
	# 				np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), 
	# 				c_int, 
	# 				c_double, 
	# 				c_double, 
	# 				np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'))
	# beta = np.array([1.5, 1.01])
	# t0 = time.time()
	# for i in range(n):
	# 	for j in range(len(probs_full)):
	# 		b = b_empty[j]
	# 		ll3c(b, 
	# 				probs, 
	# 				conc, 
	# 				len(probs), 
	# 				funmin,
	# 				1e6,
	# 				beta)
	# 		b_full[j] = b
	# t1 = time.time()
	# print(f"LL3 C: {(t1-t0)/n*1e6:.3f} us")
	
	# print(b)


	ll2ca = cloglik.ll2_array_min
	ll2ca.argtypes = (c_int, #number of probs/trial
						c_int, #number of iters
						np.ctypeslib.ndpointer(dtype=np.float64, #minima
												ndim=2, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #prob array
												ndim=2, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #conc
												ndim=1, flags='C_CONTIGUOUS'),
						np.ctypeslib.ndpointer(dtype=np.float64, #func vals
												ndim=1, flags='C_CONTIGUOUS'),
						c_double)


	niters, prob_ct = probs_full.shape
	funmin = np.zeros(niters)
	ll2ca(prob_ct,
					niters,
					b2_full, #must be size niters * 2
					probs_full, #must be size niters * prob_ct
					conc, #must be length prob_ct
					funmin,
					1e62)
	print(funmin)
	print(b2_full)

	ll3ca = cloglik.ll3_array_min
	ll3ca.argtypes = (c_int, #number of probs/trial
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
												ndim=1, flags='C_CONTIGUOUS'))


	niters, prob_ct = probs_full.shape
	funmin = np.zeros(niters)
	ll3ca(prob_ct,
					niters,
					b3_full, #must be size niters * 3
					probs_full, #must be size niters * prob_ct
					conc, #must be length prob_ct
					funmin,
					1e6,
					np.array([1.5,1.01]))
	print(funmin)
	print(b3_full)

	# print(f"LL3 C: {(t1-t0)/n*1e6:.3f} us; Min: {[b_full[i] for i in range(len(b_full))]}; Val: {funmin}")



	ll23aAIC = cloglik.array_ll2_ll3_AIC
	ll23aAIC.argtypes = (c_int, #number of probs/trial
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
												ndim=1, flags='C_CONTIGUOUS'))


	niters, prob_ct = probs_full.shape
	funmin = np.zeros(niters)
	ll23aAIC(prob_ct,
					niters,
					b3_full, #must be size niters * 3
					probs_full, #must be size niters * prob_ct
					conc, #must be length prob_ct
					funmin,
					1e6,
					np.array([1.5,1.01]))

	print(funmin)
	print(b3_full)

























	# ll3c = cfunc.ll3
	# b_type = c_double * len(b)
	# data_type = c_double * len(probs)
	# weib_type = c_double * 2
	# ll3c.restype = c_double

	# n = 200
	# t0 = time.time()
	# for i in range(n):
	# 	# ll3(np.array(b), np.array(probs), np.array(conc))
	# 	res = minimize(ll3, b, args = (probs, conc), 
	# 				method = 'BFGS', jac = ll3_jac)
	# t1 = time.time()
	# print(f"LL3 Numpy: {(t1-t0)/n*1e3:.3f} ms; Val: {res.x}")

	# def wll3c(b, probs, conc, lp):
	# 	return ll3c(b_type(*b), probs, conc, lp)

	# ll3jc = cfunc.ll3j
	# ll3jc.restype = POINTER(b_type)
	# g = b_type(*([1]*3))
	# def wll3gc(b, probs, conc, lp):
	# 	ll3jc(b_type(*b), probs, conc, lp, g)
	# 	return [g[0], g[1], g[2]]

	# print(f"LL3 C: {(t1-t0)/n*1e3:.3f} ms; Grad: {wll3gc(b, data_type(*probs), data_type(*conc), c_int(len(probs)))}")
	# print(f"LL3 C: {(t1-t0)/n*1e3:.3f} ms; Val: {wll3c(b, data_type(*probs), data_type(*conc), c_int(len(probs)))}")

	# print(f"LL3 Numpy: {(t1-t0)/n*1e6:.3f} us; Grad: {ll3_jac(np.array(b), np.array(probs), np.array(conc))}")
	# print(f"LL3 Numpy: {(t1-t0)/n*1e6:.3f} us; Val: {ll3(np.array(b), np.array(probs), np.array(conc))}")

	# t0 = time.time()
	# for i in range(n):
	# 	# ll3(np.array(b), np.array(probs), np.array(conc))
		# res = minimize(wll3c, b, args = (data_type(*probs), data_type(*conc), c_int(len(probs))), 
		# 			method = 'BFGS', jac = wll3gc)
	# t1 = time.time()
	# print(f"LL3 Numpy: {(t1-t0)/n*1e3:.3f} ms; Val: {res.x}")



	# t0 = time.time()
	# for i in range(n):
	# 	ll3c(b_type(*b), data_type(*probs), data_type(*conc), c_int(len(probs)))
	# t1 = time.time()
	# print(f"LL3 C: {(t1-t0)/n*1e6:.3f} us; Val: {ll3c(b_type(*b), data_type(*probs), data_type(*conc), c_int(len(probs)))}")




	# n = 1000
	# t0 = time.time()
	# for i in range(n):
	# 	ll3(np.array(b), np.array(probs), np.array(conc))
	# t1 = time.time()

	# t0 = time.time()
	# for i in range(n):
	# 	ll3c(b_type(*b), data_type(*probs), data_type(*conc), c_int(len(probs)))
	# t1 = time.time()
	# print(f"LL3 C: {(t1-t0)/n*1e6:.3f} us; Val: {ll3c(b_type(*b), data_type(*probs), data_type(*conc), c_int(len(probs)))}")

	# t0 = time.time()
	# b = [1.00,-0.991]
	# for i in range(n):
	# 	ll2(np.array(b), np.array(probs), np.array(conc))
	# t1 = time.time()
	# print(f"LL2 Numpy: {(t1-t0)/n*1e6:.3f} us; Val: {ll2(np.array(b), np.array(probs), np.array(conc))}")

	# ll2c = cfunc.ll2
	# ll2c.restype = c_double
	# t0 = time.time()
	# for i in range(n):
	# 	ll2c(b_type(*b), data_type(*probs), data_type(*conc), c_int(len(probs)))
	# t1 = time.time()
	# print(f"LL2 C: {(t1-t0)/n*1e6:.3f} us; Val: {ll2c(b_type(*b), data_type(*probs), data_type(*conc), c_int(len(probs)))}")

	# print(f"ll3_jac Numpy: {(t1-t0)/n*1e6:.3f} us; Val: {ll3_jac(np.array(b), np.array(probs), np.array(conc))}")

	# ll3jc = cfunc.ll3j
	# b_type = c_double * len(b)
	# data_type = c_double * len(probs)
	# ll3jc.restype = POINTER(b_type)
	# g = b_type(*([0]*3))

	# ll3jc(b_type(*b), data_type(*probs), data_type(*conc), c_int(len(probs)), g)

	# print(f"ll3_jac C: {(t1-t0)/n*1e6:.3f} us; Val: {[g[i] for i in range(len(b))]}")


	# b = [.80,-1.091]
	# print(f"ll2_jac Numpy: {(t1-t0)/n*1e6:.3f} us; Val: {ll2_jac(np.array(b), np.array(probs), np.array(conc))}")

	# ll2jc = cfunc.ll2j
	# b_type = c_double * len(b)
	# data_type = c_double * len(probs)
	# ll2jc.restype = POINTER(b_type)
	# g = b_type(*([0]*2))

	# ll2jc(b_type(*b), data_type(*probs), data_type(*conc), c_int(len(probs)), g)

	# print(f"ll2_jac C: {(t1-t0)/n*1e6:.3f} us; Val: {[g[i] for i in range(len(b))]}")

