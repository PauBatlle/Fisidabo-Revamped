import numpy as np
import matplotlib.pyplot as plt
from funcions_3 import *
import pandas as pd

#Part física - constants (SI)

g = 9.81     #obvi, no?
R = 0.04     #radi pilota
m = 0.8      #pes pilota
i = 2./5.    #factor del moment d'inercia (I/mR^2)
gamma = i/(1+i)              #factor que apareix en la força final (esfera = 2/7)


#retorna una taula amb els valors de l'acceleració predits, i també l'spin
def calcula3(taula, c, k, sigma, mu, rnu, ws_inicial):
    ws = np.zeros(len(taula)+ 1)
    ws[0] = ws_inicial
    res = []
    
    for i in range(len(taula)):
        ax,ay,ws[i+1] = ebola3r(c,k,sigma,mu,rnu,ws[i], 1/180, taula["x"][i],taula["y"][i],taula["vx"][i],taula["vy"][i])
        res.append([ax,ay])
        
    res = np.array(res)
    taula["axe"] = res[:,0]
    taula["aye"] = res[:,1]
    taula["ws"] = ws[0:-1]
    return taula




def ebola3r(c,k,sigma,mu,rnu,ws_inicial,timestep,x,y,vx,vy):
    nu = rnu/R
    sigma2 = sigma**2            #per comoditat
    #l'adam ha tingut una molt bona idea
    c = -c 
    k = -k

    #Condicions inicials (Compte!! el que sé realment a partir de les mesures són r_o i v_o, (i només en el pla x-y))
    xp0, yp0, vpx0, vpy0 = numeric(x,y,vx,vy, c, k, sigma2)
    
    #Càlculs per introduïr les cond. inicials
    r_p = [xp0, yp0, f(radius(xp0,yp0,sigma2),c,k)]
    v_p = [vpx0, vpy0, df(radius(xp0,yp0,sigma2),c,k)*vr(xp0,yp0,vpx0,vpy0,sigma2)]

    N, N_t = unit_normals(xp0,yp0,vpx0,vpy0, c, k, sigma2)
    r_o = r_p + R * N   
    v_o = v_p + R * N_t

    a_o, vws = acc_modif(r_p, r_o, v_p, v_o, ws_inicial,c,k,sigma2, mu, nu) 
    
    ws = ws_inicial + vws*timestep
        
    return a_o[0], a_o[1], ws
