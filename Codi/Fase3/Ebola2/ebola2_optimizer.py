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

def optimitza(metode):
	if metode == 1:
		d2 = minimize(nova_fitness,[c,k,sigma,mu], method='Nelder-Mead', options={'disp': True})
		print(d2)
	if metode == 2:
		bo = BayesianOptimization(lambda c, k, sigma, mu: -nova_fitness([c,k,sigma,mu]), {'c': (1.0, 3.0), 'k': (0.1, 0.8), 'sigma': (1, 2), 'mu': (0.002, 0.01)})
		gp_params = {'kernel': None,
             'alpha': 1e-5}

		# Run it again with different acquisition function
		bo.maximize(init_points = 200, n_iter=100, acq='ei', **gp_params)

		#bo.maximize(init_points=50, n_iter=25, kappa=2)
		print(bo.res['max'])
		print(bo.res['all'])


optimitza(2)