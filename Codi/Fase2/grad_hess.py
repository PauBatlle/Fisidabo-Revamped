def gradf(x,y):
    return np.array([dfx(x,y)[0], dfy(x,y)[0]]).reshape(2)

def Hessf(x,y):
    return np.array([[dfxx(x,y)[0], dfxy(x,y)[0]], [dfyx(x,y)[0], dfyy(x,y)[0]]]).reshape(2,2)