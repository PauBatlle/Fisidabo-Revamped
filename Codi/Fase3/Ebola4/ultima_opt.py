import numpy as np 
import pandas as pd 
from funcions_taula import *
from fit import *
from scipy.optimize import minimize
from numpy.random import uniform as u
from tqdm import tqdm as t
from IPython import embed
#exp3ck
c_b = (0.2,0.4)
k_b = (0.6, 0.8)
mu_b = (0.01,0.015)
s_b = (1,1.1)
def generat():
    return [u(*c_b),u(*k_b),u(*s_b),u(*mu_b)]

#x = [0.2, 0.7, 1.1,0.012]

x = generat()
y = f2exp(x, None)
i = 1
while y > 10:
	x = generat()
	y = f2exp(x, None)
	print(i, end = '\r')
	i+= 1
print(minimize(f2exp, x0 = generat(), method = "Nelder-Mead", args = {"Nfeval":0}, options={'maxiter': None, 'maxfev': None}))

"""
x"s = []
res = []
for i in t(range(5000)):
	x = generat()
	y = f3expck(x, None)
	xs.append(x)
	res.append(y)

np.save("xNITspin2", np.array(xs))
np.save("resNITspin2", np.array(res))

embed()
"""
#print(minimize(f2, x0 = x , method = "Nelder-Mead", args = {"Nfeval":0}, options={'maxiter': None, 'maxfev': None}))
