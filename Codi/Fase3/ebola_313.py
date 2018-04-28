import numpy as np
import matplotlib.pyplot as plt
from funcions_3 import *
import pandas as pd
from exporter import exporta
from distancia import distancia

#Part física - constants (SI)

g = 9.81     #obvi, no?
R = 0.04     #radi pilota
m = 0.8      #pes pilota
i = 2./5.    #factor del moment d'inercia (I/mR^2)
gamma = i/(1+i)              #factor que apareix en la força final (esfera = 2/7)

"""
#Paràmetres
c = 2     #escalat vertical en valor absolut!!!
k = 2     #escalat radial (invers)
sigma = 1.0    #elipticitat
mu = 0.01
nu = 0.5/R
ws_inicial = 0 
#Part temporal
steps = 2501
tmax = 10.   #en segons, ha de ser un float
"""

def ebola3(c,k,sigma,mu,rnu,ws_inicial,steps,timestep,x,y,vx,vy,take,escriu = False):
	nu = rnu/R
	sigma2 = sigma**2            #per comoditat
	#l'adam ha tingut una molt bona idea
	c = -c 
	k = -k

	tmax = timestep*(steps - 1)
	time = np.linspace(0, tmax, steps)

	r_p = np.zeros((steps, 3))    #posició punt contacte
	v_p = np.zeros((steps, 3)) 
	a_p = np.zeros((steps, 3)) 

	r_o = np.zeros((steps, 3))    #posició centre de masses
	v_o = np.zeros((steps, 3)) 
	a_o = np.zeros((steps, 3)) 

	ws = np.zeros(steps)          #rotació de spin de la pilota
	vws = np.zeros(steps)

	#Condicions inicials (Compte!! el que sé realment a partir de les mesures són r_o i v_o, (i només en el pla x-y))
	xp0, yp0, vpx0, vpy0 = numeric(x,y,vx,vy, c, k, sigma2)
	ws[0] = ws_inicial

	#Càlculs per introduïr les cond. inicials
	r_p[0] = [xp0, yp0, f(radius(xp0,yp0,sigma2),c,k)]
	v_p[0] = [vpx0, vpy0, df(radius(xp0,yp0,sigma2),c,k)*vr(xp0,yp0,vpx0,vpy0,sigma2)]

	N, N_t = unit_normals(xp0,yp0,vpx0,vpy0, c, k, sigma2)
	r_o[0] = r_p[0] + R * N   
	v_o[0] = v_p[0] + R * N_t

	a_o[0], vws[0], uou, aiai = acc_o(r_p[0], v_o[0], ws[0],c,k,sigma2, mu, nu)   #only per representar (nop, si que s'usa)

	#print v_p[0] - uou
	#print v_o[0] - aiai

	v_oaprox = np.zeros((steps, 3)) #UI UI

	for t in range (1, steps):
	    #Canvis de nom per a les condicions prèvies / pseudo step1
	    r_p1 = r_p[t-1] 
	    r_o1 = r_o[t-1]
	    ws1 = ws[t-1]
	    v_p1 = v_p[t-1]
	    v_o1 = v_o[t-1]
	    a_o1, vws1 = a_o[t-1], vws[t-1]
	    
	    #Steps en sí
	    r_p2, r_o2, v_p2, v_o2, a_o2, ws2, vws2 = RungeKuttastep(r_p1, r_o1, v_p1, v_o1, ws1, v_p1, v_o1, a_o1, vws1, 0.5*timestep, c, k, sigma2, mu, nu)
	        
	    r_p3, r_o3, v_p3, v_o3, a_o3, ws3, vws3 = RungeKuttastep(r_p1, r_o1, v_p1, v_o1, ws1, v_p2, v_o2, a_o2, vws2, 0.5*timestep, c, k, sigma2, mu, nu)
	    
	    r_p4, r_o4, v_p4, v_o4, a_o4, ws4, vws4 = RungeKuttastep(r_p1, r_o1, v_p1, v_o1, ws1, v_p3, v_o3, a_o3, vws3, timestep, c, k, sigma2, mu, nu)
	    
	    #Càlcul de les noves condicions
	    r_p[t] = r_p1 + (timestep/6.) * (v_p1 + 2*v_p2 + 2*v_p3 + v_p4)
	    r_o[t] = r_o1 + (timestep/6.) * (v_o1 + 2*v_o2 + 2*v_o3 + v_o4)
	    ws[t] = ws1 + (timestep/6.) * (vws1 + 2*vws2 + 2*vws3 + vws4)
	    v_o[t] = v_o1 + (timestep/6.) * (a_o1 + 2*a_o2 + 2*a_o3 + a_o4) 
	    a_o[t], vws[t], v_p[t], v_o[t] = acc_o(r_p[t], v_o[t], ws[t], c, k, sigma2, mu, nu) 

	W = np.transpose([time, ws, vws])
	N_ = np.zeros((steps, 3))
	N_t = np.zeros((steps, 3))
	for t in range (0, steps):
	    N_[t], N_t[t] = unit_normals(r_p[t,0], r_p[t,1], v_p[t,0], v_p[t,1], c, k, sigma2)

	columnes = ['time','ws','vws','rpx','rpy','rpz','vpx','vpy','vpz','x','y','roz','vox','voy','voz','aox','aoy','aoz','Nx','Ny','Nz','Ntx','Nty','Ntz']
	
	dataf = pd.DataFrame(data = np.concatenate((W, r_p, v_p, r_o, v_o, a_o, N_, N_t), 1), columns = columnes)
		
	if escriu:
		dicc = {1043:1, 1328:2, 1203:3, 649:4}
		fitness = distancia(dataf, "Dades_takes/"+str(take)+"Cart.csv")[0]
		exporta(dataf, "eb3"+"-take"+str(take)+"-fitness = "+str(fitness)+"-c="+str(c)+"-k="+str(k)+"-sigma="+str(sigma)+"-mu="+str(mu)+"-Rnu="+str(rnu)+"-ws_0="+str(ws_inicial))
	return dataf
