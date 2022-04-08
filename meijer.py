
import mpmath
import pandas as pd

def meijerg(alpha,L,mu,y):
    lista=[]
    A1 = [1]; A2 = []
    B1 = [alpha, L]; B2 = [0]
    for i in range(0,len(y)):
        z = (alpha*L*y[i])/mu[i]
        M = mpmath.meijerg([A1, A2], [B1, B2], z)
        M = round(M, 5)
        lista.append(M)
        bd = pd.DataFrame(lista)
        #print(S)
    return(bd)

