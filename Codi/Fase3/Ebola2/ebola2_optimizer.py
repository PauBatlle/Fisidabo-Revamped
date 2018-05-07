import numpy as np
from tqdm import tqdm as t
from IPython import embed
import pandas as pd
import scipy
from llegeix_escriu import esc
from fit import nova_fitness
from random import uniform
from bayes_opt import BayesianOptimization
from scipy.optimize import minimize
from random import uniform as u
#'c': (1.0, 3.0), 'k': (0.1, 0.8), 'sigma': (1, 2), 'mu': (0.002, 0.01)
c_real, k_real = 1/0.526,0.31
nova_nova_fitness = lambda x: nova_fitness([c_real,k_real,x[0],x[1]])

def optimitza(metode):
	if metode == 1:
		c1 = np.linspace(1.0,3.0,50)[0]
		k1 = np.linspace(0.1,0.8,50)[10]
		[c,k,sigma,mu] = [c1,k1,1.1,0.013]
		d2 = minimize(nova_fitness, [c,k,sigma,mu], args=(), 
			method='Nelder-Mead', tol=None, callback= print, 
			options={'disp': False, 'initial_simplex': None, 
			'maxiter': None, 'xatol': 0.0001, 'return_all': False, 
			'fatol': 0.0001, 'maxfev': None})
		#d2 = minimize(nova_fitness,[c,k,sigma,mu], method='Nelder-Mead', options={'disp': True, })
		print(d2)
	
	if metode == 2:
		bo = BayesianOptimization(lambda c, k, sigma, mu: -nova_fitness([c,k,sigma,mu]), {'c': (1.0, 3.0), 'k': (0.1, 0.8), 'sigma': (1, 2), 'mu': (0.002, 0.01)})
		# Run it again with different acquisition function
		bo.maximize(init_points = 200, n_iter=100)

		#bo.maximize(init_points=50, n_iter=25, kappa=2)
		print(bo.res['max'])
		print(bo.res['all'])

	if metode == 3:
		#Brute Force
		minim = 1e8
		M = np.zeros((30, 30))
		i = 0
		for mu in t(np.linspace(0.002,0.02,30)):
			j = 0
			for sigma in np.linspace(0,1,30):
				d = nova_fitness([1,1,sigma,mu])
				M[i][j] = min(d,99)
				j+= 1
				if d < minim:
					minim = d
					print(d, "es el nou minim")
			i+=1
		np.save("matriu_ckfixes", M)

	if metode == 4:
		[sigma,mu] = [1,0]
		d2 = minimize(nova_nova_fitness, [sigma,mu], args=(), 
			method='Nelder-Mead', tol=None, callback= print, 
			options={'disp': False, 'initial_simplex': None, 
			'maxiter': None, 'xatol': 0.0001, 'return_all': False, 
			'fatol': 0.0001, 'maxfev': None})
		#d2 = minimize(nova_fitness,[c,k,sigma,mu], method='Nelder-Mead', options={'disp': True, })
		print(d2)


optimitza(4)