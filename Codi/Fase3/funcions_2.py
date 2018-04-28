import numpy as np

g = 9.81
m = 1.

def Radius_tilla(x, y, s2):
    return np.sqrt(x**2 + s2*y*y)
    
def Modul(v):
    return np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def Normalize(v):
    return v/Modul(v)
  
def F(radi, c, k):
    return -k * np.exp(- c * radi)

def cosa(x):
    return 2*x

def funcio(x):
    return 2*x


def D_F (z, c):
    return -c * z
    
def D_D_F (z, c):
    return c * c * z

def V_r (x, v, s2):
    return (x[0]*v[0] + s2*x[1]*v[1])/Radius_tilla(x[0], x[1], s2)
    
def V_z (x, v, c, s2):
    return D_F(x[2], c) * V_r (x, v, s2)
    
def Normal(x, v, c, k, s2):
    r = Radius_tilla(x[0], x[1], s2)
    v_r = (x[0]*v[0] + s2*x[1]*v[1])/r
    dF = D_F(F(r,c,k),c)
    ddF = D_D_F(F(r,c,k),c)
    
    #càlcul i normalització del vector normal a la superfície
    n = np.array([-(dF/r) * x[0],-(dF/r) * x[1]*s2, 1])
    modul_n = Modul(n)
    n = n / modul_n          #Normalitzem el vector normal  
    
    #Component normal de l'acceleració
    a_n = ((ddF*r - dF)*v_r**2 + dF*(v[0]**2 + s2*v[1]*v[1]))/(r * modul_n) 

    #mòdul i vector de força normal
    m_normal = m*a_n + m*g/modul_n
    normal = m_normal * n
    
    return m_normal, normal
    
def Friction(v, N, mu):
    u_v = v / Modul(v) #és el vector velocitat unitari
    return - mu * N * u_v

def Acc(x, v, c, k, s2, mu):
    N, normal = Normal(x, v, c, k, s2) #mòdul i vector
    gravity = np.array([0., 0., -m*g])
    friction = Friction(v, N, mu)
    return normal + gravity + friction
    
    
def E_kinetic(v):
    v_square = v[0]**2 + v[1]**2 + v[2]**2
    return m * v_square/2.
    
def E_potential(x):
    return m * g * x[2]
    

    
    
    
