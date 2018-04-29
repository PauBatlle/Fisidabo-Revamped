import numpy as np
import pandas as pd
from funcions_2 import *
from exporter import exporta
from distancia import distancia
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

g = 9.81
m = 1.

def complete_show (tableau):
    fig = plt.figure('pos-accel-time')
    axe = fig.gca(projection='3d')
    axe.plot(tableau['x'],tableau['y'],tableau['t'], label='trajectòria')
    axe.plot(tableau['ax'],tableau['ay'],tableau['t'], label='acceleració')
    axe.plot(tableau['axe'],tableau['aye'],tableau['t'], label='accel_èbola')
    plt.xlabel('x var')
    plt.ylabel('y var')
    axe.set_zlabel('time')

    axe.legend()
    plt.show()

def distancia_acc (tableau, plot = 0):
    difx = pd.Series.as_matrix(tableau['ax']-tableau['axe'])
    dify = pd.Series.as_matrix(tableau['ay']-tableau['aye'])
    if plot:
        plt.plot(tableau['t'], difx, label = 'x')
        plt.plot(tableau['t'], dify, label = 'y')
        plt.legend()
    return (np.dot(difx,difx) + np.dot(dify,dify))/len(tableau['t'])

def calcula(taula, c, k, sigma, mu):
    k = [ebola2r(c,k,sigma,mu,taula["x"][i],taula["y"][i],taula["vx"][i],taula["vy"][i]) for i in range(len(taula))]
    k = np.array(k)
    taula["axe"] = k[:,0]
    taula["aye"] = k[:,1]
    return taula


def ebola2r(c, k, sigma, mu, x, y, vx, vy):
    #sigma:
    s2 = sigma**2
    pos = [x, y, F(Radius_tilla(x, y, s2), c, k)]
    vel = [vx, vy, 0]
    vel[2] = V_z(pos, vel,c,s2)
    acc = Acc(pos,vel,c,k,s2,mu)
    return acc