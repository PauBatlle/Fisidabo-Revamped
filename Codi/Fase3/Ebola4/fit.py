import numpy as np 
import pandas as pd 
from ebola4 import ebola4


def dq(a,b,factor = 1):
    #embed()
    a1 = a[:,0]
    a2 = a[:,1]
    b1 = b['x'].as_matrix()
    b2 = b['y'].as_matrix()
    return np.linalg.norm(a1-b1)**2 + np.linalg.norm(a2-b2)**2


dades = [pd.read_csv("../Data/"+str(i)+"Cart.csv") for i in range(1,5)]
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

#Hi ha 6 noves fitness a partir d'ara
#ebola4(ebola2, llit_exponencial, c, k, sigma, mu, nu, spininicial, xg, yg, vxg, vyg, length)

def fitness2exp(v):
    """fitness ebola2, llit exp """
    [c, k, sigma, mu] = v
    part_comuna = [True,True,c,k,sigma,mu,0.0,0.0]
    A = ebola4(*part_comuna, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    return resultat

def fitness2real(v):
    """fitness ebola2, llit real"""
    [mu] = v
    part_comuna = [True,False,-1,-1,-1,mu,0.0,0.0]
    A = ebola4(*part_comuna, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    return resultat

def fitness2expck(v):
    """fitness ebola2, llit model exponencial amb c, k reals"""
    [sigma, mu] = v
    part_comuna = [True,True,1.8,0.44,sigma,mu,0.0,0.0]
    A = ebola4(*part_comuna, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    return resultat

def fitness3exp(v):
    """fitness ebola3, llit exp """
    [c, k, sigma, mu, nu, spininicial] = v
    part_comuna = [False,True,c,k,sigma,mu,nu,spininicial]
    A = ebola4(*part_comuna, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    return resultat

def fitness3real(v):
    """fitness ebola3, llit real"""
    [mu, nu, spininicial] = v
    part_comuna = [False,False,-1,-1,-1,mu,nu,spininicial]
    A = ebola4(*part_comuna, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    return resultat

def fitness3expck(v):
    """fitness ebola3, llit model exponencial amb c, k reals"""
    [sigma, mu, nu, spininicial] = v
    part_comuna = [False,True,1.8,0.44,sigma,mu,nu,spininicial]
    A = ebola4(*part_comuna, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    return resultat