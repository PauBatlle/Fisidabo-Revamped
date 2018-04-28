import pandas as pd

def exporta(dataframe, nom):
	dataframe.to_csv('Dades_exportades/'+nom+'.csv', header = True, index = False)
