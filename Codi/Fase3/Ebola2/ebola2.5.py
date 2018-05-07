import numpy as np
import pandas as pd
from funcions_2 import *
g = 9.81
m = 1.

"""
steps = 1500 
timestep = 1/180

mu = 0.005
c = 2.2
k = 0.4
x = 0
y = 0.5
vx = 0.8
vy = 0
s = 1.5
"""

def ebola2(c, k, sigma, mu, steps, timestep, x, y, vx, vy, take, escriu = False):
    #sigma:
    s2 = sigma**2
    pos = np.zeros((steps, 3)) 
    vel = np.zeros((steps, 3)) 
    acc = np.zeros((steps, 3))
    energy = np.zeros((steps))
    error = np.zeros((steps))
    temps_maxim = timestep*(steps-1)
    time = np.linspace(0, temps_maxim, steps)
    pos[0] = [x, y, F(Radius_tilla(x, y, s2), c, k)]
    vel[0] = [vx, vy, 0]
    vel[0][2] = V_z(pos[0], vel[0], c, s2)
    acc[0] = Acc(pos[0],vel[0], c,k,s2,mu)
    for t in range (1, steps):
        #RUNGE-KUTTA
        x1 = pos[t-1]
        v1 = vel[t-1]
        a1 = Acc(x1,v1,c,k,s2,mu)
        x2 = x1 + 0.5*v1*timestep
        v2 = v1 + 0.5*a1*timestep
        a2 = Acc(x2,v2,c,k,s2, mu)
        x3 = x1 + 0.5*v2*timestep
        v3 = v1 + 0.5*a2*timestep
        a3 = Acc(x3,v3,c,k,s2,mu)
        x4 = x1 + v3*timestep
        v4 = v1 + a3*timestep
        a4 = Acc(x4,v4,c,k,s2,mu)        
        pos[t] = x1 + (timestep/6.)*(v1+2*v2+2*v3+v4)
        vel[t] = v1 + (timestep/6.)*(a1+2*a2+2*a3+a4)
        acc[t] = Acc(pos[t],vel[t],c,k,s2,mu) #only per representar
        
    for t in range (0, steps):
        x = pos[t]
        energy[t] = E_kinetic(vel[t]) + E_potential(pos[t])
        error[t] = x[2] - F(Radius_tilla(x[0], x[1], s2), c, k)

    o = np.column_stack((time, pos[:,0], pos[:,1], pos[:,2], vel[:,0], vel[:,1], vel[:,2], acc[:,0], acc[:,1], acc[:,2], energy, error))
    columnes = ['time','x','y','z','vx','vy','vz','ax','ay','az', 'energy', 'error']
    dataf = pd.DataFrame(data = o, columns = columnes)
    return dataf
