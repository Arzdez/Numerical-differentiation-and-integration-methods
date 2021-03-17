import numpy as np
import matplotlib.pyplot as plt

#sins = [ np.sin(i) for i in np.arange(0,np.pi/2,1)]
#x = [i for i in np.arange(0,np.pi/2,1)]

def sins(x):
    return np.sin(x)
a = 0
b = np.pi/2

# правые прямоугольники
def right(a,b,n):
    sumy = 0
    dx=(b-a)/n
    result = 0
    for i in range(n):
        sumy+=dx
        result += sins(sumy)

    return(result*dx)

def midle(a,b,n):
    dx = (a-b)/n
    xd = dx/2
    result = 0
    sumy = 0
    sumy2 = 0
    for i in range(n):
        sumy += dx
        sumy2 = sumy - xd
        result += sins(sumy2)
    return result*dx

def trapeciya(a,b,n):
    dx = (a-b)/n
    result = 0
    sumy = 0
    xd = 0
    for i in range(n):
        xd+=dx
        sumy+=sins(xd)
    
    I = dx*( ((sins(a)+sins(b))/2)+sumy)
    return( I)

def simpson(a,b,n):
    dx = (a-b)/n
    xd = dx/2
    sumy = 0
    sumy2 = 0
    suma = 0
    for i in range(1,n*2-1):
        sumy += dx
        sumy2 = sumy - xd
        suma +=sins(sumy2)
    
    #print(len(T))
    #sum1 = 4*sum([sins(T[i]) for i in np.arange(1,n*2+1,1)])
    #sum2 = 2*sum([sins(T[i]) for i in np.arange(1,n*2,1)])

    #print(sum1)
    #print(sum2)
    #I = dx/3 * ( sins(a)+sins(b)+sum1+sum2)
   # return I

simpson(a,b,100)

#print(right(a,b,100))
#print(midle(a,b,100))
#print(trapeciya(a,b,100))
#clprint(simpson(a,b,100))