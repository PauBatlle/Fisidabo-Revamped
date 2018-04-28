#LES FUNCIONS QUE UTILITZA l'EBOLA 3.1.2
import numpy as np
from numpy.linalg import norm

#Part física - constants (SI)

g = 9.81     #obvi, no?
R = 0.04     #radi pilota
m = 0.8      #pes pilota
i = 2./5.    #factor del moment d'inercia (I/mR^2)
gamma = i/(1+i)              #factor que apareix en la força final (esfera = 2/7)


def radius(x,y,sigma2):
    return np.sqrt(x*x + sigma2*y*y)

#Velocitat radial: necessita la posició x y i les seves derivades, retorna la velocitat radial
def vr(x, y, vx, vy, sigma2):
    return np.dot([x, y*sigma2],[vx, vy])/radius(x,y,sigma2)

#f FUNCIÓ ALTURA (RADIAL)
def f(r,c,k):
    return c * np.exp(k*r)
def df(r,c,k):
    return k*c * np.exp(k*r)
def ddf(r,c,k):
    return k*k*c * np.exp(k*r)

#la funció F(x,y) és f(r(x,y))
#DERIVADES PARCIALS DE F
def derivades_F(x, y, c, k, sigma2):
    r = radius(x,y, sigma2)
    grad_r = np.array([x, sigma2*y])/r
    hess_r = sigma2 * np.array([[y*y, - x*y],
                       [- y*x, x*x]])/r**3
    
    df = k * f(r,c,k) #derivade de f respecte r
    ddf = k * df  #2a derivade
    
    grad_F = df * grad_r
    hess_F = ddf * np.transpose([grad_r]) * grad_r + df * hess_r
    return grad_F, hess_F


#Vector unitari normal i la seva derivada temporal
def unit_normals(x, y, vx, vy, c, k, sigma2):
    grad_F, hess_F = derivades_F(x, y, c, k, sigma2)
    n_ = np.append(-grad_F, 1)      #vector normal a la superfície
    n = np.sqrt(np.dot(n_, n_))     #el seu mòdul
    N = n_/n
    n_t = np.array([-np.dot(hess_F[0], [vx, vy]), 
                    -np.dot(hess_F[1], [vx, vy]),
                   0])              #derivada del vector normal respecte el temps

    N_t = n_t/n - N/n * np.dot(n_t, N)
    
    return N, N_t

#Troba la part ortogonal a n del vector v
def ortogonalitza(v, n):
    return v - np.dot(v,n)*n

#Calcula les forces sobre el centre de masses O, en retorna l'acceleració
#També fa tangent la velocitat del c.m. a la superfície, i calcula la velocitat del desplaçament de 
# p (punt on es produeix el contacte de la bola amb el terra) en l'espai. 
# (alerta, p no és el punt de la bola en contacte amb el terra, que està instantàniament quiet)

def acc_o(r_p, v_o, ws, c, k, sigma2, mu, nu): #els que porten barra baixa són vectors, els altres, nombres  (per a un cert moment t)
    #Calculem coses bàsiques
    
        #PART no temporal
    (x,y,f) = r_p
    grad_F, hess_F = derivades_F(x, y, c, k, sigma2)
    n_ = np.append(-grad_F, 1)       #vector normal a la superfície
    n = np.sqrt(np.dot(n_, n_))      #el seu mòdul
    N_ = n_/n                        #vector normal unitari         #OBS: vectors en majúscula són unitaris,
                                                                    #vectors amb " _ " tenen 3 components (com ja dit)
    
        #PART 1a derivada temporal
    v_o = ortogonalitza(v_o, N_)     #Correcció ortogonalitat de v_o
    v_p = vel_p(grad_F, hess_F, v_o) #Càlcul de la velocitat en la superfície del punt de contacte(no del punt de l'esfera/sòlid)
    #Alternativament, calcula la derivada discreta, usant: N_t2  aprox=  N2-N1
    #Haurem de demanar d'entrada el N_ previ (podria entrar r_pi)
    
    (vx,vy,vf) = v_p
    n_t = np.array([-np.dot(hess_F[0], [vx, vy]), 
                    -np.dot(hess_F[1], [vx, vy]),
                   0])               #derivada del vector normal respecte el temps
    vo = np.sqrt(np.dot(v_o,v_o))
    V = v_o/vo
    
    #CALCULEM FORCES
        #Part Normal
    an = - np.dot(v_o, n_t)/n    #acceleració normal, del centre de masses
    Fn = m*(g*N_[2] + an)            #força normal, feta pel llit elàstic (apunt, mgN_[2] és el negatiu de la projecció de 
                                    # la força de la gravetat sobre la normal)  
                                    # (no cal calcular-la)
                                    # (potser millor calcular-la a banda, i fer-ne plot)
           
        #Part Tangencial
    g_t = -g * np.array([0,0,1] - N_*N_[2])           #acceleració tangencial causada per la gravetat (N_[2] = N*[0,0,1] = 1/n)
    F_t = gamma*m * (R*ws * np.cross(N_, n_t)/n - g_t) #Força "aparent" deguda a la rotació de la pilota,
                                                                      #i la seva fricció amb el terra
    a_t = F_t/m + g_t - mu*Fn*V  #acceleració tangencial, del centre de masses, potser era més ràpid  computacionalment 
                                    #el mètode anterior
    
    #Obtenim acceleració
    a_o = a_t + an * N_
    
    #Calculem canvi spin
    vws = np.dot(n_t, np.cross(N_, v_o))/(R*n)
    
    if (abs(ws) > 1e-3):
         vws = vws - nu*ws/abs(ws)
            
    return a_o, vws, v_p, v_o

#Retorna v_p en funció de r_p i v_o 
# (velocitat del punt de contacte, en funció de la seva posició i la velocitat del centre de masses)
def vel_p(grad_F, hess_F, v_o):
    n = np.sqrt(np.dot(grad_F, grad_F) + 1)
    a, b, c = [-1 - grad_F[1]**2, grad_F[0]* grad_F[1], -1 - grad_F[0]**2]
    
    #solucionem un sistema d'equacions: v_p + R * N_t = v_o (incògnites vx, vy, sudem de la tercera component) 
    #1a equació A*vx + B*vy = v_o[0]
    #2a equació C*vx + D*vy = v_o[1]
    [[A, B], [C,D]] = (R/n**3) * np.dot([[a,b], [b,c]], hess_F) + np.identity(2)
    
    #Solucionem per gauss: (no seidel) (alerta C i B peetits)
    vy = (v_o[1] - C/A * v_o[0] )/(D - C*B/A)
    vx = (v_o[0] - B*vy)/A 
    
    vf = np.dot(grad_F, [vx, vy]) #¡¡¡¡¡Aquest pas marca el motiu de perquè v_p es comporta tant bé en l'eix z!!!!!
                                    #podria calcular-ho per gauss també la tercera component i veure
    return np.array([vx, vy, vf])


def comprova(r_p, v_p, r_g, v_g, c, k, sigma2): #Retorna cert si la solució de numèric és l'esperada
    r_p = [r_p[0], r_p[1], f(radius(r_p[0],r_p[1], sigma2),c,k)]
    v_p = [v_p[0], v_p[1], df(radius(r_p[0],r_p[1], sigma2),c,k)*vr(r_p[0],r_p[1],v_p[0],v_p[1],sigma2)]
    N, N_t = unit_normals(r_p[0],r_p[1],v_p[0],v_p[1], c, k, sigma2)
    nr_g = r_p + R*N
    nv_g = v_p + R*N_t
    return (max(norm(nr_g[0:2]-r_g), norm(nv_g[0:2]-v_g)) < 1e-10)



def numeric(xg0,yg0,vgx0,vgy0, c, k, sigma2):
    #Càlculs per introduïr les cond. inicials
    r_g = [xg0, yg0]
    v_g = [vgx0, vgy0]
    N, N_t = unit_normals(xg0,yg0,vgx0,vgy0, c, k, sigma2)
    r_p = r_g - R * N[0:2]
    v_p = v_g - R * N_t[0:2]
    norma_dif = 1
    niter = 0
    while norma_dif > 1e-10 and niter < 30: #a_r guarda la de la iteració anterior
        niter += 1
        ar_p = r_p
        av_p = v_p
        N, N_t = unit_normals(ar_p[0],ar_p[1],av_p[0],av_p[1], c, k, sigma2)
        r_p = r_g - R * N[0:2]
        v_p = v_g - R * N_t[0:2]
        norma_dif = max(norm(r_p-ar_p), norm(v_p-av_p))

    if niter == 30: 
        raise Exception("El mètode no ha convergit")

    if not comprova(r_p, v_p, r_g, v_g, c, k, sigma2):
        raise Exception("El mètode no ha convergit a un bon lloc")

    return r_p[0], r_p[1], v_p[0], v_p[1]



#ªªªªªªªªªªªªªªª#
# Runge - Kutta #
#ªªªªªªªªªªªªªªª#

#Fa un pas de Runge Kutta (se'n fan 4 o 3 vegades més que passos habituals, aquesta funció en concret, es crida 3 vegades)
#i significa inicial, ac significa actual
#Obs: v_pi no s'usa, però es passa per notació/continuïtat... com es vulgui dir, coherència
def RungeKuttastep(r_pi, r_oi, v_pi, v_oi, wsi, v_pac, v_oac, a_oac, vwsac, timestep, c, k, sigma2, mu, nu):
    r_o2 = r_oi + v_oac*timestep
    r_p2 = r_pi + v_pac*timestep
    ws2 = wsi + vwsac*timestep
    v_o2 = v_oi + a_oac*timestep
    a_o2, vws2, v_p2, v_o2 = acc_o(r_p2, v_o2, ws2, c, k, sigma2, mu, nu)
    
    return r_p2, r_o2, v_p2, v_o2, a_o2, ws2, vws2   #El 2 vol dir que és el de sortida


#########################################
#Funció acceleració per a l'èbola retall#
#########################################
def acc_modif(r_p, r_o, v_p, v_o, ws, c, k, sigma2, mu, nu): #els que porten barra baixa són vectors, els altres, nombres  (per a un cert moment t)
    #Calculem coses bàsiques
    
        #PART no temporal
    (x,y,f) = r_p
    grad_F, hess_F = derivades_F(x, y, c, k, sigma2)
    n_ = np.append(-grad_F, 1)       #vector normal a la superfície
    n = np.sqrt(np.dot(n_, n_))      #el seu mòdul
    N_ = n_/n                        #vector normal unitari         #OBS: vectors en majúscula són unitaris,
                                                                    #vectors amb " _ " tenen 3 components (com ja dit)
    
    (vx,vy,vf) = v_p
    n_t = np.array([-np.dot(hess_F[0], [vx, vy]), 
                    -np.dot(hess_F[1], [vx, vy]),
                   0])               #derivada del vector normal respecte el temps
    vo = np.sqrt(np.dot(v_o,v_o))
    V = v_o/vo
    
    #CALCULEM FORCES
        #Part Normal
    an = - np.dot(v_o, n_t)/n    #acceleració normal, del centre de masses
    Fn = m*(g*N_[2] + an)            #força normal, feta pel llit elàstic (apunt, mgN_[2] és el negatiu de la projecció de 
                                    # la força de la gravetat sobre la normal)  
                                    # (no cal calcular-la)
                                    # (potser millor calcular-la a banda, i fer-ne plot)
           
        #Part Tangencial
    g_t = -g * np.array([0,0,1] - N_*N_[2])           #acceleració tangencial causada per la gravetat (N_[2] = N*[0,0,1] = 1/n)
    F_t = gamma*m * (R*ws * np.cross(N_, n_t)/n - g_t) #Força "aparent" deguda a la rotació de la pilota,
                                                                      #i la seva fricció amb el terra
    a_t = F_t/m + g_t - mu*Fn*V  #acceleració tangencial, del centre de masses, potser era més ràpid  computacionalment 
                                    #el mètode anterior
    
    #Obtenim acceleració
    a_o = a_t + an * N_
    
    #Calculem canvi spin
    vws = np.dot(n_t, np.cross(N_, v_o))/(R*n)
    
    if (abs(ws) > 1e-3):
         vws = vws - nu*ws/abs(ws)
            
    return a_o, vws
