import numpy as np 
import pandas as pd 
from ebola4 import ebola4
from scipy.optimize import minimize

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
'''
mu = 0,inf
k = 0.1,0.8
sigma = 1,2
mu = 0,inf
nu = 0, inf
ws = -5,5
'''

ROP = 250 #mes que aixo pot voler dir rop
#FUNCIONS PER A LA OPTIMITZACIO

def analitzarop(num, traj):
    moduls = [np.linalg.norm([traj[i,0],traj[i,1]]) for i in range(len(traj))]
    print(num, "maxim a prop: ", min(moduls))

def f2exp(v, info, ana = False): #4 PARAMS: (C,K,SIGMA,MU)
    """fitness ebola2, llit exp """
    [c, k, sigma, mu] = v
    part_comuna = [True,True,c,k,sigma,mu,0.0,0.0]
    A = ebola4(*part_comuna, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    
    if info is not None:    

        if info['Nfeval'] == 0:
            print ('{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s}  {5:9s}'.format('Iter', 'c', 'k', 'sigma', 'mu', 'f(X)'))

        if info['Nfeval']%5 == 0:
            print ('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}  {4:3.6f}   {5:3.6f}'.format(info['Nfeval'], c, k, sigma, mu, resultat))
        info['Nfeval'] += 1

    if ana:
        if resultat > ROP or np.isnan(resultat):
            if np.isnan(resultat):
                print("NAN!")
            print(v, resultat)
            print(dq(A, d1p),dq(B,d2p),dq(C,d3p),dq(D,d4p))
            analitzarop(1,A)
            analitzarop(2,B)
            analitzarop(3,C)
            analitzarop(4,D)
            print("-------")

    return resultat

def f2real(v, info, ana = False): #1 PARAM: (MU)
    """fitness ebola2, llit real"""
    [mu] = v
    part_comuna = [True,False,-1,-1,-1,mu,0.0,0.0]
    A = ebola4(*part_comuna, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    
    if info is not None:
        if info['Nfeval'] == 0:
            print ('{0:4s}   {1:9s}   {2:9s}'.format('Iter', 'mu', 'f(X)'))

        if info['Nfeval']%5 == 0:
            print ('{0:4d}   {1: 3.6f}   {2: 3.6f}'.format(info['Nfeval'], mu, resultat))
        info['Nfeval'] += 1
    
    if ana:
        if resultat > ROP or np.isnan(resultat):
            if np.isnan(resultat):
                print("NAN!")
            print(v, resultat)
            print(dq(A, d1p),dq(B,d2p),dq(C,d3p),dq(D,d4p))
            analitzarop(1,A)
            analitzarop(2,B)
            analitzarop(3,C)
            analitzarop(4,D)
            print("-------")

    return resultat

def f3exp(v, info, ana = False): #9 PARAMS: (C,K,SIGMA,MU,NU,SPININICIAL*4)
    """fitness ebola3, llit exp """
    [c, k, sigma, mu, nu, s1, s2, s3, s4] = v
    spininicial = [s1,s2,s3,s4]
    part_comuna = [False,True,c,k,sigma,mu,nu]
    A = ebola4(*part_comuna, s1, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, s2, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, s3, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, s4, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))

    if info is not None:
        # display information
        if info['Nfeval'] == 0:
            print ('{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s}  {5:9s} {6:9s}   {7:9s}  {8:9s}  {9:9s}  {10:9s}'.format('Iter', 'c', 'k', 'sigma', 'mu', 'nu', 's1', 's2', 's3', 's4', 'f(X)'))

        if info['Nfeval']%5 == 0:
            print ('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}  {4:3.6f}    {5:3.6f}    {6:3.6f}   {7:3.6f} {8:3.6f} {9:3.6f} {10:3.6f}'.format(info['Nfeval'], c, k, sigma, mu, nu, s1, s2, s3, s4, resultat))
        info['Nfeval'] += 1

    if ana:
        if resultat > ROP or np.isnan(resultat):
            if np.isnan(resultat):
                print("NAN!")
            print(v, resultat)
            print(dq(A, d1p),dq(B,d2p),dq(C,d3p),dq(D,d4p))
            analitzarop(1,A)
            analitzarop(2,B)
            analitzarop(3,C)
            analitzarop(4,D)
            print("-------")

    return resultat

def f3real(v, info, ana = False): #6 PARAMS: (MU, NU, SPININICIAL*4)
    """fitness ebola3, llit real"""
    [mu, nu, s1, s2, s3, s4] = v
    part_comuna = [False,False,-1,-1,-1,mu,nu]
    A = ebola4(*part_comuna, s1, x01, y01, vx01, vy01, len(d1p)-1)
    B = ebola4(*part_comuna, s2, x02, y02, vx02, vy02, len(d2p)-1)
    C = ebola4(*part_comuna, s3, x03, y03, vx03, vy03, len(d3p)-1)
    D = ebola4(*part_comuna, s4, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(A, d1p)+dq(B,d2p)+dq(C,d3p)+dq(D,d4p))
    
    if info is not None:
        # display information
        if info['Nfeval'] == 0:
            print ('{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s}'.format('Iter', 'mu', 'nu', 's0', 'f(X)'))

        if info['Nfeval']%5 == 0:
            print ('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}  {4:3.6f}'.format(info['Nfeval'], mu, nu, spininicial, resultat))
        info['Nfeval'] += 1
    
    if ana:
        if resultat > ROP or np.isnan(resultat):
            if np.isnan(resultat):
                print("NAN!")
            print(v, resultat)
            print(dq(A, d1p),dq(B,d2p),dq(C,d3p),dq(D,d4p))
            analitzarop(1,A)
            analitzarop(2,B)
            analitzarop(3,C)
            analitzarop(4,D)
            print("-------")


    return resultat


#FUNCIONS PER A LA OPTIMITZACIO AMB C = 1.8 FIXA

def f2expc(v, info, ana = False): #3 PARAMS: (K, SIGMA, MU)
    """fitness ebola2, llit exp, c fixa """
    return f2exp(np.array([1.8]+list(v)), info, ana)
def f2expck(v, info, ana = False): #2 PARAMS: (SIGMA, MU)
    """fitness ebola2, llit exp, c fixa """
    return f2exp(np.array([1.8, 0.44]+list(v)), info, ana)

def f3expck(v, info, ana = False): #7 PARAMS: (SIGMA, MU, NU, SPININICIAL*4)
    """fitness ebola3, llit exp, c fixa """
    return f3exp(np.array([1.8,0.44]+list(v)), info, ana)
