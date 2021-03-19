import numpy as np
import matplotlib.pyplot as plt
import time

startTime = time.time()
#sins = [ np.sin(i) for i in np.arange(0,np.pi/2,1)]
#x = [i for i in np.arange(0,np.pi/2,1)]

def sins(x):
    return np.sin(x)
a = 0
b = np.pi/2
# правые прямоугольники
def Iright(Fx,a,b,n):
    #sumy = 0
    dx=(b-a)/n
    result = 0
    for i in range(n):
        result += Fx(dx*i)

    return(result*dx)

def Imidle(Fx,a,b,n):
    dx = (b-a)/n
    xd = dx/2
    result = 0
    sumy = 0
    sumy2 = 0
    for i in range(n):
        sumy += dx
        sumy2 = sumy - xd
        result += Fx(sumy2)
    return result*dx

def Itrapeciya(Fx,a,b,n):
    dx = (b-a)/n
    sumy = 0
    for i in range(n):
        sumy+=Fx(dx*i)
    I = dx*( ((Fx(a)+Fx(b))/2)+sumy)
    return( I)

def Isimpson(Fx,a,b,n): 
    h=((b-a)/n)
    s=Fx(a)+Fx(b)
    for i in range(1,n):
        if i%2==0:
            s+= 2*Fx(a+h*i) 
        else:
            s+= 4*Fx(a+h*i)
    return(s*h/3)

#def ImonteKarlo(Fx,a,b,n):
    #p1
    

print(Iright(sins,a,b,1000))
print(Imidle(sins,a,b,1000))
print(Itrapeciya(sins,a,b,1000))
print(Isimpson(sins,a,b,1000))


endTime = time.time()
totalTime = endTime - startTime
print("Время = ", totalTime)