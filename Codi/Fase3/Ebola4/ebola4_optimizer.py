import numpy as np
from tqdm import tqdm as t
from fit import *
from funcions_taula import *
from random import uniform
from bayes_opt import BayesianOptimization
from scipy.optimize import minimize, brute
from IPython import embed
#funcions = [fitness2exp,fitness2real,fitness2expck, fitness3exp, fitness3real, fitness3expck]
#condicions_in = [[1.8,0.44,1,0.001],[0.001],[1,0.001],[1.8,0.44,1,0.001,0.0,0.0],[0.001,0,0],[1,0.001,0,0]]

mu_b = (0, 1)
#print(minimize(f2expc, x0 = [0.4,1.1,0.05], args = ({'Nfeval':0}), bounds = bounds, options={'ftol': 0.0001, 'maxiter': INF}))


v1 = np.linspace(*mu_b, 100)
#v2 = np.linspace(*nu_b, 15)
#v3 = np.linspace(*spin_b, 15)

M = np.zeros(100)
for i in t(range(100)):
	M[i] = f2real([v1[i]], None)

embed()