import pandas as pd

def distancia(set1,set2):
    vector = [(set1["x"][i]-set2["x"][i])**2 + (set1["y"][i]-set2["y"][i])**2 for i in range(min(len(set1), len(set2)))]
    s = sum(vector)
    l = len(vector)
    return s, l, s/l
