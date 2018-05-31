import numpy as np
from tqdm import tqdm as t
from fit import *
from random import uniform
from bayes_opt import BayesianOptimization
from scipy.optimize import minimize, brute
from IPython import embed
#funcions = [fitness2exp,fitness2real,fitness2expck, fitness3exp, fitness3real, fitness3expck]
#condicions_in = [[1.8,0.44,1,0.001],[0.001],[1,0.001],[1.8,0.44,1,0.001,0.0,0.0],[0.001,0,0],[1,0.001,0,0]]

#print(minimize(f2expc, x0 = [0.1,1.1,0.05], args = ({'Nfeval':0})))

sigma_b = (1, 1.5)
mu_b = (0, 1)
#print(minimize(f2expc, x0 = [0.4,1.1,0.05], args = ({'Nfeval':0}), bounds = bounds, options={'ftol': 0.0001, 'maxiter': INF}))

#v1 = np.linspace(*k_b, 10)
v1 = np.linspace(*sigma_b, 20)
v2 = np.linspace(*mu_b, 20)

M = np.zeros((20,20))
for i in t(range(20)):
	for j in t(range(20)):
			M[i,j] = f2expck([v1[i],v2[j]], None)

embed()