import numpy as np
from tqdm import tqdm as t
from fit import *
from funcions_taula import *
from random import uniform
from bayes_opt import BayesianOptimization
from scipy.optimize import *
from IPython import embed
#funcions = [fitness2exp,fitness2real,fitness2expck, fitness3exp, fitness3real, fitness3expck]
#condicions_in = [[1.8,0.44,1,0.001],[0.001],[1,0.001],[1.8,0.44,1,0.001,0.0,0.0],[0.001,0,0],[1,0.001,0,0]]

#print(minimize(f2expc, x0 = [0.1,1.1,0.05], args = ({'Nfeval':0})))
"""
INF = int(1e8)
c_b = (1.2,2.4)
k_b = (0.1,0.8)
sigma_b = (1,1.5)
mu_b = (0, 0.2)
bounds = (k_b,sigma_b,mu_b)

print(minimize(f2expc, x0 = [0.177777,1.05556,0.022222], method = "Powell", args = ({"Nfeval":0})))
"""
sigma_b = (1,1.5)
mu_b = (0, 0.5)
nu_b = (0, 0.5)
spin_b = (-2, 2)

v0 = np.linspace(*sigma_b, 5)
v1 = np.linspace(*mu_b, 5)
v2 = np.linspace(*nu_b, 5)
v3 = np.linspace(*spin_b, 5)

M = np.zeros((5,5,5,5))

for l in t(range(5)):
	for i in t(range(5)):
		for j in t(range(5)):
			for k in t(range(5)):
				M[l,i,j,k] = f3expck([v0[l],v1[i],v2[j],v3[k]], None)

embed()