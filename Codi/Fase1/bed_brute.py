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
n = 4

M = np.zeros((10,10))
v1 = np.linspace(50, 400, len(M))
v2 = np.linspace(0.5, 4, len(M[0]))
for i in t([2*n,2*n+1]):
	for j in t(range(len(M[0]))):
		M[i][j] = fitness([v1[i],v2[j]])

np.save("llits_grid"+str(n), M)