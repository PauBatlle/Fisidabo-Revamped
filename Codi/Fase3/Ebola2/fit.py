import numpy as np
import sympy as s
import pandas as pd
import matplotlib.pyplot as plt
from ebola_2 import ebola2
import distancia
from importlib import reload
from distancia import distancia as d
from llegeix_escriu import esc
from IPython import embed
dades = [pd.read_csv("../Data/"+str(i)+"Cart.csv") for i in range(1,5)]

def dq(a,b,factor = 1):
    #embed()
    a1 = a[::factor]['x'].as_matrix()
    a2 = a[::factor]['y'].as_matrix()
    b1 = b['x'].as_matrix()
    b2 = b['y'].as_matrix()
    return np.linalg.norm(a1-b1)**2 + np.linalg.norm(a2-b2)**2

d1 = dades[0]
x01 = d1["x"][175]
y01 = d1["y"][175]
vx01 = d1["vx"][175]
vy01 = d1["vy"][175]
d2 = dades[1]
x02 = d2["x"][270]
y02 = d2["y"][270]
vx02 = d2["vx"][270]
vy02 = d2["vy"][270]
d3 = dades[2]
x03 = d3["x"][100]
y03 = d3["y"][100]
vx03 = d3["vx"][100]
vy03 = d3["vy"][100]
d4 = dades[3]
x04 = d4["x"][100]
y04 = d4["y"][100]
vx04 = d4["vx"][100]
vy04 = d4["vy"][100]
d1p = dades[0].iloc[175:838]
d2p = dades[1].iloc[270:1158]
d3p = dades[2].iloc[100:1063]
d4p = dades[3].iloc[100:429]

def nova_fitness(v, factor = 1):
    [c,k,sigma,mu] = v
    ebolaa = ebola2(c, k, sigma, mu, len(d1p)*factor, (1/180)/factor,x01,y01,vx01,vy01,"dummy",False)[['x','y']]
    ebolab = ebola2(c, k, sigma, mu, len(d2p)*factor, (1/180)/factor,x02,y02,vx02,vy02,"dummy",False)[['x','y']]
    ebolac = ebola2(c, k, sigma, mu, len(d3p)*factor, (1/180)/factor,x03,y03,vx03,vy03,"dummy",False)[['x','y']]
    ebolad = ebola2(c, k, sigma, mu, len(d4p)*factor, (1/180)/factor,x04,y04,vx04,vy04,"dummy",False)[['x','y']]
    A = dq(ebolaa, d1p, factor)
    B = dq(ebolab,d2p, factor)
    C = dq(ebolac,d3p, factor)
    D = dq(ebolad,d4p, factor)
    resultat = 0.25*(A+B+C+D)
    #esc("results.csv", [c, k, sigma, mu, resultat])
    #print(resultat, A, B, C, D)
    return resultat