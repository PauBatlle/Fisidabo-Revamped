import numpy as np
from tqdm import tqdm as t
from fit import *
from random import uniform
from bayes_opt import BayesianOptimization
from scipy.optimize import minimize, brute
from IPython import embed
import argparse
"""
parser = argparse.ArgumentParser()
parser.add_argument("i")
args = parser.parse_args()
"""
#print(minimize(f2expc, x0 = [0.1,1.1,0.05], args = ({'Nfeval':0})))

v = [[1.02631578947, 0.157894736842],
[1.05263157895, 0.157894736842], 
[1.05263157895, 0.263157894737], 
[1.05263157895, 0.473684210526], 
[1.07894736842, 0.210526315789], 
[1.07894736842, 0.421052631579], 
[1.13157894737, 0.210526315789], 
[1.13157894737, 0.315789473684], 
[1.18421052632, 0.210526315789], 
[1.21052631579, 0.210526315789], 
[1.23684210526, 0.263157894737], 
[1.31578947368, 0.263157894737], 
[1.34210526316, 0.473684210526], 
[1.36842105263, 0.421052631579], 
[1.39473684211, 0.368421052632], 
[1.39473684211, 0.421052631579], 
[1.44736842105, 0.210526315789], 
[1.44736842105, 0.368421052632], 
[1.44736842105, 0.473684210526], 
[1.47368421053, 0.210526315789], 
[1.47368421053, 0.263157894737]]

"""
if args.i is not None:
	i = int(args.i)
else:
	i = 0
print(i)
"""


print(minimize(f3expck, x0 = [1,0.005,0.005,-5] , method = "Nelder-Mead", args = {"Nfeval":0}, options={'maxiter': None, 'maxfev': None}))


#scipy.optimize.minimize
#(fun, x0, args=(), method='Nelder-Mead', 
#tol=None, callback=None, 
#options={'func': None, 'maxiter': None, 
#'maxfev': None, 'disp' : False, 'return_all': False, 
#'initial_simplex': None, 'xatol': 0.0001, 
#'fatol': 0.0001, 'adaptive': False})
"""
for i in range(len(v)):
	print(i)
	print("---------------------")
	print(minimize(f2expck, x0 = [v[i]], method = "Nelder-Mead", args = None, options={'maxiter': None, 'maxfev': None}))
	print("---------------------")
#print(minimize(f2expck, x0 = [1.0,0.0], method = "Nelder-Mead", args = ({'Nfeval':0}), options={'maxiter': None, 'maxfev': None}))
#print(minimize(f3expc, x0 = [ks[i],sigmes[i],mus[i],nus[i],s0[i]], method = "Nelder-Mead", args = ({'Nfeval':0}), options={'maxiter': None, 'maxfev': None}))
"""
"""
rranges = (slice(1, 2, 0.1), slice(0, 0.2, 0.02))
from scipy import optimize
s = optimize.brute(f2expck, rranges, full_output=True)
print(s)
embed()
"""