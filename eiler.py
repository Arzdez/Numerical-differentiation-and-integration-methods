import math
import numpy as np
import matplotlib.pyplot as plt
import random
def F(y):
    return y

def G(x,y,A,W,t):
    return -y-x**3+A*np.sin(W*t)

def eiler(F,G,N,A,w,h):
    X,Y = np.empty(N+1), np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()

    for i in range(0,N):

        X[i+1] = X[i]+h*F(Y[i])

        Y[i+1] = Y[i]+h*G(X[i],Y[i],A,w,h*(i))
    
    return(X,Y)

D, B =  eiler(F,G,100000,420,3.9,0.0001)

plt.plot(D[:],B[:])
plt.show()