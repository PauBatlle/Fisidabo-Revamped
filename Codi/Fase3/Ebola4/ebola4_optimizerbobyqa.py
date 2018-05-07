import numpy as np
from tqdm import tqdm as t
from fit import *
import pybobyqa

#funcions = [fitness2exp,fitness2real,fitness2expck, fitness3exp, fitness3real, fitness3expck]
#condicions_in = [[1.8,0.44,1,0.001],[0.001],[1,0.001],[1.8,0.44,1,0.001,0.0,0.0],[0.001,0,0],[1,0.001,0,0]]
"""
mu = 0.002 - 0.02
c = 1.0 - 3.0
k = 0.1 - 0.8
s = 1 - 2
mu = [0.002, 0.02]
nu = [0.1, 0.8]/R
ws_inicial = [-5, 5] 
"""
lower = np.array([0.002])
upper = np.array([0.025])
funcions = [fitness2real]
condicions_in = np.array([[0.002]])
for i in range(len(funcions)):
	#print(minimize(funcions[i], condicions_in[i], method='Nelder-Mead', options={'maxiter': None, 'maxfev': None, 'disp': True}))
	print(pybobyqa.solve(funcions[i], condicions_in[i], bounds = (lower,upper 