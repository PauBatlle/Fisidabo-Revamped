import pandas as pd
from os import listdir as ls

def llegeix_dades():
    dades = {}
    for arxiu in ls("Dades_takes"):
    	t = pd.read_csv('Dades_takes/'+arxiu)
    	dades[arxiu[0:2]] = t[["t","x","y","vx","vy","ax","ay"]]
    return dades