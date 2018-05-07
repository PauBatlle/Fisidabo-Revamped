import numpy as np
from tqdm import tqdm as t
from fit import *
from random import uniform
from bayes_opt import BayesianOptimization
from scipy.optimize import minimize, brute

#funcions = [fitness2exp,fitness2real,fitness2expck, fitness3exp, fitness3real, fitness3expck]
#condicions_in = [[1.8,0.44,1,0.001],[0.001],[1,0.001],[1.8,0.44,1,0.001,0.0,0.0],[0.001,0,0],[1,0.001,0,0]]

funcions = [fitness2real]
condicions_in = [0.2]
"""
for i in range(len(funcions)):
	print(minimize(funcions[i], condicions_in[i], method='Nelder-Mead', options={'maxiter': None, 'maxfev': None, 'disp': True}))
"""
brute_force = True
if brute_force:
	rranges = (slice(0.001,0.03,(0.03-0.001)/10), slice(0,0.1,0.01),slice(-0.05,0.05,0.01))
	print(brute(fitness3real, rranges))