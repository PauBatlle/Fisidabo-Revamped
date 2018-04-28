
# Simulació llit elàstic
# Aquest codi es basa en discretitzar un llit elàstic similar al del fisidabo, que és una xarxa quadrada, i aplicar en ell les equacions del moviment. Cada punt de la xarxa quadrada se li associa una petita massa 'dm'. I la corda es considera com una sèrie de molles molt rígides, amb una distància de repòs 'dist_repos', i una constant de rigidesa 'K'. Molles entre els punts d'unió de la xarxa que estan units per corda.
# 
# Per a aclariments mirar una de les llibretes que regalaven de la UPC i que té altres coses variades de mates.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm as t
from IPython import embed
from scipy.interpolate import interp1d
import pandas as pd
import scipy
import os
import time
from llegeix_escriu import esc
#Valors inicials: 175, 1.4

def calcula_un_llit(elasticitat, tensio_inicial, log = False, save = False):

	#print("He entrat amb parametres", elasticitat, tensio_inicial)
	"""
	if os.path.exists("Function-Results/"+string):
		return np.load("Function-Results/"+string+"/tensor.npy"), np.load("Function-Results/"+string+"/X.npy"), np.load("Function-Results/"+string+"/Z.npy") 
	"""
	string = "EL="+str(elasticitat)+"-TI="+str(tensio_inicial)

	if True:
		Elasticitat_cordill = elasticitat
		Tensio_inicial = tensio_inicial
		Elasticitat_molla_metall = 400 #Paràmetre que en un principi no cal modificar


		# ### El control de la sortida en un .txt es controla al final del programa

		#PROPIETATS DEL TEIXIT/CORIDLL
		x_nodes = 214 #Nombre de creuaments entre cordes quan ens movem en l'eix x, que és el llarg
		y_nodes = 150 #Nombre de creuaments entre cordes quan ens movem en l'eix y, que és el curt

		#Mida del llit elàstic sense gravetat (sense deformacio, tal i com ho fa en començar la simulació), en cada dimensió.
		#En metres. Només de la part que és corda, no les molles que subjecten
		x_mida_inici = 3.6 
		y_mida_inici = 1.8

		#Considerem de moment una elasticitat i una massa uniforme:
		densitat = 0.01 #En Kg per metre de corda
		K = Elasticitat_cordill #Mòdul d'elasticitat de la corda, multiplicat per la secció. En Newtons (per metre estirat i en un metre de corda).


		#PROPIETATS DE LES MOLLES DE METALL
		K_molla = Elasticitat_molla_metall #Elasticitat de les molles de metall al voltant del llit (són les que permeten més deformació i lliguen el llit 
		# a les barres de suport laterals)
		L_molla = 0.40 #Longitud molles de metall
		#Podríem afegir també la densitat de les molles


		#SITUACIO INICIAL
		Tensio_inicial = Tensio_inicial #En Newtons, força que hi ha sobre el llit elàstic en ambdós eixos en començar el llit horizontal.

		#Penjem una anella de massa de les següents característiques:
		massa_esfera = 67.5 #En kg, massa penjada
		Radi = 0.1
		Radi_intern = 0.05

		#Domini temporal
		#timestep = 0.0005
		#temps_simulació = 1 #En segons
		g = 9.8 #Vity


		#REINICIAT DEL LLIT
		#A continuació elegim les posicions inicials dels creuaments i dels punts importants (fixació als voltants/contorns), 
		#i ho posem en la matriu Pos
		x_indexs = x_nodes+4
		y_indexs = y_nodes+4
		Pos = np.zeros((x_indexs, y_indexs, 3))   #Les posicions de cada un dels creuaments en l'espai de 3 dimensions, 
		#començant per 2 i acabant per x_nodes+1 o y_nodes+1 depenent de la direcció. En la posició 0 i en la posició final x_nodes+3 
		#o y_nodes+3 correspon als punts de l'espai que són fixes per a aquest llit elàstic, els que s'enganxen a les barres 
		#metàl·liques laterals en el nostre cas

		#Nombre de masses que tindran dinàmica en cada cordill segons l'eix
		x_masses = x_nodes+2
		y_masses = y_nodes+2
		Gird = np.mgrid[0:x_mida_inici:(x_masses)*1j,0:y_mida_inici:(y_masses)*1j]  #crea una graella, com les que serveixen 
		#per a fer gràfics de superfícies, amb la part complexa de divisions i la part real d'interval a dividir.
		#Aquesta Gird seria la part de teixit del llit elàstic

		Pos[1:-1,1:-1,0] = Gird[0]  #Assignem les posicions en l'eix X de la part de texit
		Pos[1:-1,1:-1,1] = Gird[1]  #Assignem les posicions en l'eix Y de la part de texit
		#En l'eix Z la posició comença amb un 0

		#Falta afegir la longitud de les molles de metall.
		Pos[1:-1,1:-1,0:2] = Pos[1:-1,1:-1,0:2] + L_molla

		#I omplir els espais de les posicions fixes de contorn que havíem deixat buits
		Pos[0,:,1] = Pos[1,:,1]
		Pos[-1,:,1] = Pos[1,:,1]
		Pos[:,-1,1] = Pos[:,-2,1] + L_molla
		Pos[:,0,0] = Pos[:,1,0]
		Pos[:,-1,0] = Pos[:,1,0]
		Pos[-1,:,0] = Pos[-2,:,0] + L_molla

		#Mides en "repòs" (sense gravetat), hi pot haver tensió en repòs no nul·la
		longitud_inici_x = Pos[2,1,0] - Pos[1,1,0]
		longitud_inici_y = Pos[1,2,1] - Pos[1,1,1]

		#CINEMÀTICA INICIAL
		Vel = np.zeros((x_masses,y_masses,3)) #Velocitat inicialitzada a 0

		altura_centre = np.array([]) #Inicialitzem un vector que guardarà l'altura del centre de llit elàstic

		# Plot the surface.
		#PLOT A SURFACE
		"""
		fig = plt.figure('Llit')
		ax = fig.gca(projection='3d')

		fig.suptitle("Llit elàstic")
		plt.xlabel("X (meters)")
		plt.ylabel("Y (meters)")
		ax.set_zlabel('Height (meters)')
		ax.set_zlim(-1,0.1)

		surf = ax.plot_surface(Pos[:,:,0], Pos[:,:,1], Pos[:,:,2], antialiased=True)
		"""
		#plt.savefig("Llit estàtic Buit K20, 3D")

		#PREPARACIÓ DE LES FORCES DE TENSIÓ

		#Equació que ens dóna la tensió T = (K/x_eq)(long_rep - x_eq) = K(long_rep/x_eq - 1)
		# T/K +1 = Long_rep/x_eq --> x_eq = long_rep/(T/K +1)
		x_dist_repos = longitud_inici_x/(Tensio_inicial/K + 1)
		y_dist_repos = longitud_inici_y/(Tensio_inicial/K + 1)
		molla_dist_repos = L_molla/(Tensio_inicial/K_molla + 1)

		k_x = K/x_dist_repos #Constant d'elasticitat de les molles en l'eix x
		k_y = K/y_dist_repos #Constant d'elasticitat de les molles en l'eix y
		k_molla = K/molla_dist_repos #Constant d'elasticitat de les molles de metall

		deformacio_inicial_x = longitud_inici_x - x_dist_repos
		deformacio_inicial_y = longitud_inici_y - y_dist_repos
		deformacio_inicial_molla = L_molla - molla_dist_repos

		#DINÀMICA
		#Matriu de masses
		dm = densitat*(x_dist_repos + y_dist_repos) #Perquè hi ha quatre mitjos de massa de corda per cada creuament, i perquè 
		#considerem que en repòs, les densitats de les cordes en ambdues direccions són iguals i valen 'densitat'
		#Considerem ara la matriu de les masses de cada punt
		Masses = dm * np.ones((x_masses, y_masses)) #Les caselles (i,j) es refereixen a la massa del node en la posició (i+1,j+1).
		#Allí on s'enganxen les molles no es considera exactament un node (almenys no un node com a creuament de teixit).
		#L'índex comença a partir de 0 i va fins x_nodes+1 o fins y_nodes+1.

		#-------------------------------------------------------------------#
		# Afegirem una massa uniforme en forma de circumferència de 67,5 kg #
		#-------------------------------------------------------------------#
		index_centre = (int(x_masses/2), int(y_masses/2)) #Aquest índex fa referència només allí on ocorre dinàmica
		index_pos_centre = (index_centre[0]+1, index_centre[1]+1) 

		#Primer contem el nombre de creuaments en una circumferència de radi Radi centrada en index_centre
		def distancia_centre(index):
		    "Retorna la distància euclídia de l'encreuament a la que fa referència l'índex respecte l'encreuament central del llit"
		    metrica = np.array([longitud_inici_x, longitud_inici_y]) #Distancia que avancem per cada moviment de l'índex
		    index_relatiu = index - np.array (index_centre)
		    posicio_relativa = index_relatiu*metrica
		    return np.linalg.norm(posicio_relativa)


		Anella = np.zeros((x_masses, y_masses))
		for x in range(x_masses):
		    for y in range(y_masses):
		        index = np.array([x,y])
		        if (Radi_intern < distancia_centre(index) < Radi):
		            #print('dins l'anella')
		            Anella[x,y] = 1
		            
		caselles_afectades = sum(sum(Anella))
		densitat_uniforme_bola = massa_esfera/caselles_afectades
		Massa_bola = Anella * densitat_uniforme_bola
		Masses = Masses + Massa_bola

		mida_inicial_x = longitud_inici_x
		mida_inicial_y = longitud_inici_y
		
		if log:
		#A saco de data
			print("Pos centre:", index_pos_centre)
			print ("Massa d'un encreuament =", dm)
			print ('Massa del cordill =', sum(sum(Masses)) - massa_esfera) #Massa del llit elàstic
			estimació_mida_errors = 2e-16*4*K*(x_nodes+1)/(x_mida_inici*dm)
			print('Estimació_error_numèric_per_iteració:', estimació_mida_errors)
			print('Mida inicial de les molles, eixx', mida_inicial_x)
			print('Mida inicial de les molles, eixy', mida_inicial_y)
			print('Allargament inicial eix X', deformacio_inicial_x)
			print('Mida en repòs eix X', x_dist_repos)
			print('Mida en repòs eix Y', y_dist_repos)
			print('Constant elàstica eix X', K/x_dist_repos)
			print('Constant elàstica eix Y', K/y_dist_repos)
			print('Constant elàstica molla', K_molla/molla_dist_repos)



		#Una "Iteració"
		def acceleracio (Pos, Masses):
		    """Introduïm en Pos les posicions de cada creuament i lloc d'enganxament, en Masses les masses en els creuaments,
		     a partir de les altres variables globals, retorna l'acceleració que correspon a cada massa degut a les forces
		     que li fan les molles"""
		    #Genera matriu de posicions relatives:
		    x_PosRelatives = Pos[1:,1:-1,:] - Pos[:-1,1:-1,:] #Ens dóna els vectors que uneixen els diferents creuaments en l'eix x. 
		    #S'obté en fer la diferència entre les posicions de cada un d'ells i el d'una posició desplaçada cap al negatiu en l'eix x.
		    y_PosRelatives = Pos[1:-1,1:,:] - Pos[1:-1,:-1,:] #Ens dóna els vectors que uneixen els diferents creuaments en l'eix y. 
		    #S'obté en fer la diferència entre les posicions de cada un d'ells i el d'una posició desplaçada cap al negatiu en l'eix y.
		    #En la direcció en la que no fem la resta retallem una posició per sobre i una per sota, que són les que estem fent servir 
		    #per marcar les condicions de contorn, de manera que no les considerem quina distància tenen entre elles,
		    #perquè sabem que no cambiarà.

		    #Construïm la matriu de tensions en magnitud:
		        #Primer calculem les longituds de les molles (norma de la diferència de posicions -> distància):
		    x_Distancies = np.linalg.norm(x_PosRelatives, axis = 2)
		    y_Distancies = np.linalg.norm(y_PosRelatives, axis = 2)
		    
		    #A partir d'aquí hem de tractar diferent les molles de metall
		        #Mirem quina deformació han tingut
		    x_Deformacio_teixit = x_Distancies[1:-1, :] - x_dist_repos
		    x_Deformacio_molla = x_Distancies[(0, -1), :] - molla_dist_repos
		    y_Deformacio_teixit = y_Distancies[:, 1:-1] - y_dist_repos
		    y_Deformacio_molla = y_Distancies[:,(0, -1)] - molla_dist_repos
		        #Calculem tensió corresponent
		    x_Tensio_teixit = x_Deformacio_teixit * (K/x_dist_repos)
		    x_Tensio_molla = x_Deformacio_molla * (K_molla/molla_dist_repos)
		    y_Tensio_teixit = y_Deformacio_teixit * (K/y_dist_repos)
		    y_Tensio_molla = y_Deformacio_molla * (K_molla/molla_dist_repos)
		        #Ho ajuntem en unes mateixes matrius de tensions
		    x_Tensio = np.concatenate(([x_Tensio_molla[0]], x_Tensio_teixit, [x_Tensio_molla[1]]))
		    y_Tensio = np.concatenate((np.transpose([y_Tensio_molla[:,0]]), y_Tensio_teixit, np.transpose([y_Tensio_molla[:,1]])), axis = 1)


		    #print(x_Distancies-Pos[1,0,0])
		    #print(sum(x_Distancies))
		    #print(x_Tensio, np.shape(x_Tensio))
		    #print(y_Tensio)
		    
		    #A partir d'aquí ja tot és independent de què tinguem al llit, només tractem amb les forçes trobades
		    #Matriu de Tensions com a vectors força
		    x_DireccionsRelatives = x_PosRelatives/x_Distancies[:,:,np.newaxis]
		    y_DireccionsRelatives = y_PosRelatives/y_Distancies[:,:,np.newaxis]
		    
		    #print (np.linalg.norm(y_DireccionsRelatives, axis = 2))
		    
		    x_VectorTensio = x_Tensio[:,:,np.newaxis]*x_DireccionsRelatives #El vull fer aquí és multiplicar cada vector de posició relativa 
		    y_VectorTensio = y_Tensio[:,:,np.newaxis]*y_DireccionsRelatives #(paral·lel a la corresponent molla) pel coeficient de tensió que li toqui.
		        
		    #Calculem forçes resultants degudes a les cordes o pesos en cada eix.
		    x_ForçaResultant = x_VectorTensio[1:] - x_VectorTensio[:-1]
		    #print ('Forçax\n', x_ForçaResultant)
		    y_ForçaResultant = y_VectorTensio[:,1:] - y_VectorTensio[:,:-1]
		    #print ('Forçay\n', y_ForçaResultant)
		    Gravetat = - g

		    ForçaResultant = x_ForçaResultant + y_ForçaResultant #Força Resultant només del llit elàstic
		    #print ('Força\n', ForçaResultant)
		    accelResultant = ForçaResultant/Masses[:,:,np.newaxis] #Point-wise division: força sobre cada punt dividit per la massa de cada punt
		    accelResultant[:,:,2] = accelResultant[:,:,2] + Gravetat
		    #print ('Acceleració\n', accelResultant)
		    return accelResultant

		def EulerStep(Posicio_creuament, Velocitat, Acc, timestep, atenuacio = 0):
		    """Fa un step de l'algorisme d'Euler per a aproximar el moviment. Possibilitat d'afegir coeficient d'atenuacio.
		    No calcula l'acceleració"""
		    Resultat_Vel = (1-timestep*atenuacio)*Velocitat + Acc*timestep #l'atenuació causa un efecte de dampening, de fricció i amortiguament de la velocitat
		    Resultat_Pos_creuament = Posicio_creuament + Resultat_Vel*timestep

		    return Resultat_Pos_creuament, Resultat_Vel


		def RungeKuttaStep(Pos, Vel_1, timestep, atenuacio = 0):
		    Acc_1 = acceleracio(Pos, Masses)
		    Pos_copia = np.copy(Pos) #Em guardo una còpia de la posició inicial
		    Pos_creuament = Pos_copia[1:-1,1:-1] #Quan modifico aquesta porció també modifico Pos_copia? SÍ, si faig servir Pos_creuament[índex],
		    #NO si ho faig fent una igualació de tot. I quan modifico Pos_copia[1:-1,1:-1], modifico igualment aquesta porció. 
		    #Això és el que faig servir a continuació
		    [Pos_copia[1:-1,1:-1], Vel_2] = EulerStep(Pos_creuament, Vel_1, Acc_1, 0.5*timestep, atenuacio = atenuacio)
		    Acc_2 = acceleracio(Pos_copia, Masses)
		    
		    [Pos_copia[1:-1,1:-1], Vel_3] = EulerStep(Pos_creuament, Vel_2, Acc_2, 0.5*timestep, atenuacio = atenuacio)
		    Acc_3 = acceleracio(Pos_copia, Masses)
		    
		    [Pos_copia[1:-1,1:-1], Vel_4] = EulerStep(Pos_creuament, Vel_3, Acc_3, timestep, atenuacio = atenuacio)
		    Acc_4 = acceleracio(Pos_copia, Masses)
		    
		    Pos_copia[1:-1, 1:-1] = Pos[1:-1, 1:-1] + (timestep/6) * (Vel_1 + 2*Vel_2 + 2*Vel_3 + Vel_4)
		    Vel = Vel_1 + (timestep/6) * (Acc_1 + 2*Acc_2 + 2*Acc_3 + Acc_4)
		    return Pos_copia, Vel

		factor = 2
		timestep = 0.00012/factor

		#COMPROVEM SI HI HA ERRORS D'INTEGRACIÓ QUE FARAN DIVERGIR EL CODI# 
		#Alerta, per als valors petits de tensió inicial que hi ha ara, no funciona bé#

		DesplEuler = timestep*timestep #En un principi era el desplaçament en la primera iteració, però ara és un desplaçament de mostra
		#Fet amb 1 m/s^2 d'acceleració cap a un sentit

		Allargament = np.sqrt(DesplEuler**2 + mida_inicial_x**2)-mida_inicial_x
		k = K/x_dist_repos
		Tensio_inicial = k * deformacio_inicial_x
		Tensio_Euler = Tensio_inicial + k*Allargament
		Acceleració_reactiva_primerstep_per_mollax = DesplEuler*Tensio_Euler/(mida_inicial_x*dm) #La projecció de la tensió sobre l'eix z
		Acceleració_reactiva_primerstep_per_mollay = DesplEuler*Tensio_Euler/(mida_inicial_y*dm)
		#Com que els cantons estan tirats per dues molles, i suposem que al primer pas és on es produeixen les pitjors
		#acceleracions, que és mentida, perquè al segon step serà pitjor:
		Estimació_acceleració_reactiva = Acceleració_reactiva_primerstep_per_mollax + Acceleració_reactiva_primerstep_per_mollay

		if log:
			print ('Tensió inicial:', Tensio_inicial, '\nTensió després del primer Eulerstep', Tensio_Euler)
			print(DesplEuler*g, "\nFactor d'oscil·lació:", 4*Estimació_acceleració_reactiva)

			if (Estimació_acceleració_reactiva < 1/4 ):
			    #Amb que fos menor que g, de manera que no hi hagués cap punt que agafi velocitat cap amunt, en el segon step, 
			    #la cosa semblava controlada, però no, es veu que la fita que es necessita és 4 vegades menys 
			    print ('El codi no divergirà, molt probablement, per empirisme, equivocat quan les K són molt més altes que la Tensió inicial')
			else:
			    print ('ALERTA! Probablement el codi divergirà')
			    factorreducció = np.sqrt(Estimació_acceleració_reactiva*4)
			    print ('Redueix el timestep dividint-lo per un factor de:', factorreducció)

		for i in range(60*factor):
		    for j in range (200):
		        [Pos[1:-1,1:-1], Vel] = EulerStep(Pos[1:-1,1:-1], Vel, acceleracio(Pos, Masses), timestep, atenuacio = 10)
		    #print ('Resultatis euler\n',[Pos, Vel])
		    Pos_centre = Pos[index_centre]
		    altura_centre = np.append(altura_centre, Pos_centre[2])


		Pos[:,:,0] = Pos[:,:,0] - x_mida_inici/2 - L_molla
		Pos[:,:,1] = Pos[:,:,1] - y_mida_inici/2 - L_molla

	if save:
		os.makedirs("Function-Results/"+string)
		np.save("Function-Results/"+string+"/tensor", Pos)
		np.save("Function-Results/"+string+"/X", Pos[:,index_pos_centre[1],0])
		np.save("Function-Results/"+string+"/Z", Pos[:,index_pos_centre[1],2])
	return Pos, Pos[:,index_pos_centre[1],0], Pos[:,index_pos_centre[1],2]

#Posant-li True em torna el tall amb l'eix

def llegeix_reals():
	df = pd.read_csv("dadesllitreals.csv", header = None)
	df = df/100
	df = df[[4,5]]
	df = df[abs(df[4]) > 0.1]
	return df

def fitness(v):
	df = llegeix_reals()
	[elasticitat, tensio_inicial] = v
	df[5] = df[5] + 0.05
	#string = "EL="+str(elasticitat)+"-TI="+str(tensio_inicial)
	tensor_total, x_to_inter, y_to_inter = calcula_un_llit(elasticitat, tensio_inicial)
	f = interp1d(x_to_inter, y_to_inter, kind='cubic') #Chill, no es una funcio cubica sino splines cubics
	vec = [float(f(j)) for j in df[4]]
	#resultat = np.linalg.norm(vec-df[5])
	g = lambda h: np.linalg.norm(vec-(df[5]+h))
	#altura va de -0.2 a 0.2
	x0 = 0
	bounds = [[-0.5, 0.5]]
	k = scipy.optimize.minimize(g, x0 = x0, jac = False, bounds = bounds)	
	if not k.success:
		print("Bug ", v)
		print(k)
	resultat = k.fun
	altura = k.x[0]
	esc("results.csv", [elasticitat, tensio_inicial, altura, resultat])
	#print("fitness de",elasticitat, tensio_inicial,"=",resultat)
	"""
	if not os.path.isfile("Function-Results/"+string+"/fitness3.txt"):
		fitxer = open("Function-Results/"+string+"/fitness3.txt",'w')
		fitxer.write(str(resultat))
		fitxer.close()
	"""
	#print([elasticitat, tensio_inicial, altura, resultat])
	return resultat


#v_el = np.linspace(20,200,40)#[i*8:(i+1)*8]
#v_tens = np.linspace(0,8,40)