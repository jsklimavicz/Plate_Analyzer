import numpy as np
import math
from ctypes import *
import time 
from scipy.optimize import minimize, least_squares


def ll3(b, probs, conc, sigma_squared = 1e6, weibull_param=[2,1]):
	'''
	Log-likelihood function of the three parameter dose-response curve 
						    b2
				y = ------------------
					1 + exp(b0 + b1*x)
	wherein priors are b0, b1 ~ MVN(0, sigma*I2) and 
	b2 ~ Weibull(weibull_params).

	Paramerization of Weibull(x|l, k) = k/l * (x/l)^(k-1) * exp(-(x/l)**k)
	I forget why I implemented this with Weibull 
	instead of Beta, but it works fine so I haven't
	changed it. 
	'''
	b0, b1, b2 = b
	if (b2 <= 1e-10 or b2 > 1. ): return(1e10)
	xi = np.exp(b0+b1*conc)
	alpha = 1+xi # >1.0
	l = np.log(alpha)
	if (min(alpha)-b2 <= 1e-10): return(1e10)
	wk = weibull_param[0]
	wl = weibull_param[1]*(wk/(wk-1))**(1/wk)
	#note: from benchmarking, b0*b0 is approximately 3x faster than b0**2 for a float.
	ll = -(b0*b0 + b1*b1)/(2*sigma_squared) + sum(probs*np.log(alpha-b2)) +\
			math.log(b2)*sum((1-probs)) - sum(l)
	ll += -((b2/wl)**wk) + (wk-1)*math.log(b2) #+ np.log(wk) - wk*np.log(wl)

	ll2 = -(b0*b0 + b1*b1)/(2*sigma_squared)  -((b2/wl)**wk) + (wk-1)*math.log(b2)
	print(f"Before start: ll = {ll2}")

	# for i in range(len(conc)):
	# 	ll2 += probs[i]*math.log(alpha[i]-b2) + math.log(b2)*(1-probs[i]) - l[i]
	# 	# print(f"After {i}: ll = {ll2}")
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

def ll3_jac(b, probs, conc,sigma_squared = 1e6, weibull_param=[2,1]):
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
	xi = np.exp(b0+b1*conc)
	alpha = 1+xi
	d = (alpha - b2)
	m = probs*xi / d
	l = xi/alpha
	wk = weibull_param[0]
	wl = weibull_param[1]*(wk/(wk-1))**(1/wk)
	g0 = -b0/(sigma_squared) + sum(m) - sum(l)
	g1 = -b1/(sigma_squared) + sum(conc*m) - sum(conc*l)
	g2 = (wk-1)/b2 - (wk/b2)*((b2/wl)**wk) - sum(probs/d) + \
			sum((1-probs)/b2)
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

	so_file = "./ll.so"
	cfunc = CDLL(so_file)

	conc = np.array([8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.0625, 0.03125])
	probs = np.array([0, 0, 0, 0, 0.6875, 0.846153846153846, 1, 0.857142857142857, 0.923076923076923, 0.954545454545455, 0.875, 0, 0, 0, 0, 0, 0, 0, 0.083333333333333, 0.15, 0.916666666666667])
	b = [.90,-1.191,.92]

	t0 = time.time()
	t1 = time.time()
	n=1

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

	# t0 = time.time()
	# for i in range(n):
	# 	# ll3(np.array(b), np.array(probs), np.array(conc))
	# 	res = minimize(wll3c, b, args = (data_type(*probs), data_type(*conc), c_int(len(probs))), 
	# 				method = 'BFGS', jac = wll3gc)
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
	print(f"LL3 Numpy: {(t1-t0)/n*1e6:.3f} us; Val: {ll3(np.array(b), np.array(probs), np.array(conc))}")

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

