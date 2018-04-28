import numpy as np
from tqdm import tqdm as t
from IPython import embed
import pandas as pd
import scipy
from llegeix_escriu import esc
from llitelastic import *
from random import uniform

""" de llit elastic hereda la funcio fitness(v)
que va de R^2 --> R 
		(el, t0) --> R


ESPAI DE PARAMETRES SOBRE EL QUE VOLEM OPTIMITZAR:
K: 50 - 400
TensióInicial: 0.5 - 4
"""
def optimitza(metode):
	if metode == 1:
		""" Mètode 1: Scipy.minimize nelder-mead + no constriants """ 
		x0 = [175, 1.41]
		return scipy.optimize.minimize(fitness, x0, method='Nelder-Mead', tol=None, callback=None, options={'disp': True, 'initial_simplex': None, 'maxiter': None, 'xatol': 0.0001, 'return_all': False, 'fatol': 0.0001, 'maxfev': None})
	if metode == 2:
		x0 = [uniform(50, 400), uniform(0.5,4)]
		bounds = [[50, 400], [0.5,4]]
		return scipy.optimize.minimize(fitness, x0)

	if metode == 3:
		""" Mètode 3: Algorismes genètics """ 

	if metode == 4:
		""" Mètode 4: Bayesian Optimization """ 
		from bayes_opt import BayesianOptimization
		# Lets find the maximum of a simple quadratic function of two variables
		# We create the bayes_opt object and pass the function to be maximized
		# together with the parameters names and their bounds.
		bo = BayesianOptimization(lambda el, t0: -fitness([el,t0]), {'el': (50, 400), 't0': (0.5, 4)})
		# One of the things we can do with this object is pass points
		# which we want the algorithm to probe. A dictionary with the
		# parameters names and a list of values to include in the search
		# must be given.
		# Once we are satisfied with the initialization conditions
		# we let the algorithm do its magic by calling the maximize()
		# method.
		#bo.explore({'el': [112.31372046097945], 't0': [3.2815978953503064]})
		bo.maximize(init_points=10, n_iter=15, kappa=2)

		# The output values can be accessed with self.res
		# If we are not satisfied with the current results we can pickup from
		# where we left, maybe pass some more exploration points to the algorithm
		# change any parameters we may choose, and the let it run again.
		# Making changes to the gaussian process can impact the algorithm
		# dramatically.
		gp_params = {'kernel': None,
		             'alpha': 1e-5}

		# Run it again with different acquisition function
		bo.maximize(n_iter=5, acq='ei', **gp_params)

		# Finally, we take a look at the final results.
		print(bo.res['max'])
		print(bo.res['all'])	


g = optimitza(4)
print(g)
from IPython import embed
embed()