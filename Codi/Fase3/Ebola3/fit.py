import numpy as np 
import pandas as pd 
from ebola3final import ebola3
from llegeix_escriu import esc

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


def nova_fitness(v, factor = 1):
    [mu,nu,spin0] = v
    ebolaa = ebola3(mu, nu, spin0, x01, y01, vx01, vy01, len(d1p)-1)
    ebolab = ebola3(mu, nu, spin0, x02, y02, vx02, vy02, len(d2p)-1)
    ebolac = ebola3(mu, nu, spin0, x03, y03, vx03, vy03, len(d3p))
    ebolad = ebola3(mu, nu, spin0, x04, y04, vx04, vy04, len(d4p)-1)
    resultat = 0.25*(dq(ebolaa, d1p, factor)+dq(ebolab,d2p, factor)+dq(ebolac,d3p, factor)+dq(ebolad,d4p, factor))
    esc("results.csv", [mu,nu,spin0,resultat])
    return resultat
print(nova_fitness([0.4,0.2,0]))