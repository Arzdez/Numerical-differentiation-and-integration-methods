import math
import numpy as np
import matplotlib.pyplot as plt

def F(x,t):
    return 



def eiler(Fx,y0,n,N):
    h = N/n
    t =  np.linspace(0)
    y = np.zeros(N+1)
    y.append(y0)
    for i in range(1,N):
        y[i+1]=y[i]+h*Fx(y[i],h)
    return y,t
