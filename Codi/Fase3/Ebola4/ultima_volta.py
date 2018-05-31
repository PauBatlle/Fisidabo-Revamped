import numpy as np 
import pandas as pd 
from funcions_taula import *
from scipy.optimize import minimize
from numpy.random import uniform as u


c_b = (1.4, 2.2)
k_b = (0.2, 0.5)
sigma_b = (1,1.5)
mu_b = (0,0.15)
nu_b = (0,0.5)
s0_b = (-5,5)

def generat():
    return [u(*c_b),u(*k_b),u(*sigma_b),u(*mu_b),u(*nu_b),u(*s0_b)]

"""
funcions = [f1,f2,f3,f4]
x = generat()
y = f4(x, None)
i = 0
while np.isnan(y) or y > 5:
    x = generat()
    y = f4(x, None)
    print(i, end = '\r')
    i+= 1
print("--")
"""
x4 = [ 2.08050412,  0.42154183,  1.06881135,  0.00395291, -0.02453657,
       -3.3707517 ]
print(minimize(f4, x0 = x4 , method = "Nelder-Mead", args = {"Nfeval":0}, options={'maxiter': None, 'maxfev': None}))
