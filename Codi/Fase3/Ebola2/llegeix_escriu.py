import pandas as pd 

def esc(fitxer, vector_resultats):
	""" Funció per guardar resultats de fitness en un .csv """
	df = pd.read_csv(fitxer)
	df.loc[len(df)] = vector_resultats
	df.to_csv(fitxer, index = False)
