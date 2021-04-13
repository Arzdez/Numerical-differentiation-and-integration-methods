import math
import numpy as np
import matplotlib.pyplot as plt
import random
from math import log
def Ex(X):
    return 1/X
#Осцилятор уеды
def F_eiler(y):
    return y

def G_Eiler(x,y,A,W,t):
    return -y-x**3+A*np.sin(W*t)

#Осцилятор Рёсслера
def Xaxs(y,z):
    return (- y - z)

def Yaxs(x,y,a):
    return ( x +a*y)

def Zaxs(x,z,b,r):
    return ( b+(x-r)*z)

#Генератор Кияшко-Пиковского-Рабиновича
def Ff(z):
    return ( (8.592*z)-(22*z*z) + (14.408*(z*z*z)) )
def Xax(x,y,z,u,g):
    return 2*u*x + y - g*z

def Yax(x):
    return -x
def Zax(x,z):
    b =(x - Ff(z) ) / 0.2
    return b


#Решение системы уеды
def eiler_ueda(F,G,N,A,w,h):
    X,Y = np.empty(N+1), np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()

    for i in range(0,N):

        X[i+1] = X[i]+h*F(Y[i])

        Y[i+1] = Y[i]+h*G(X[i],Y[i],A,w,h*(i))
    
    return(X,Y)

def Runge_kutte_for_Ueda(F,G,N,A,w,h):
    
    X,Y = np.empty(N+1), np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()
    
    for i in range(N):
        K1X = F(Y[i])
        K1Y = G(X[i],Y[i],A,w,h*(i))

        K2X = F(Y[i]+(h/2)*K1X)
        K2Y = G(X[i]+(h/2)*K1X,Y[i]+(h/2)*K1Y,A,w,(h*i)+h/2)

        K3X = F(Y[i]+(h/2)*K2X)
        K3Y = G(X[i]+(h/2)*K2X,Y[i]+(h/2)*K2Y,A,w,(h*i)+h/2)

        K4X = F(Y[i]+(h/2)*K3X)
        K4Y = G(X[i]+(h/2)*K3X,Y[i]+(h/2)*K3Y,A,w,(h*i)+h/2)

        X[i+1] = X[i]+(h/6)*(K1X+(2*K2X)+(2*K3X)+K4X)
        Y[i+1] = Y[i]+(h/6)*(K1Y+(2*K2Y)+(2*K3Y)+K4Y)

    return X,Y

# Решение системы Рёсслера
def eiler_Ressler(N,h,a,r,b=0.2):
    X,Y,Z = np.empty(N+1), np.empty(N+1),np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()
    Z[0] = np.random.rand()

    for i in range(N):
        X[i+1] = X[i]+ h * ( -Y[i] - Z[i] )

        Y[i+1] = Y[i]+ h* ( X[i] + a * Y[i] )

        Z[i+1] = Z[i] + h * (  b + (X[i] - r)*Z[i] )
    
    return X,Y,Z

def Runge_kutte_for_ressler(Xas,Yas,Zas,N,h,a,r,b=0.2):
    X,Y,Z= np.empty(N+1), np.empty(N+1),np.empty(N+1)

    X[0] = np.random.rand()
    Y[0] = np.random.rand()
    Z[0] = np.random.rand()

    for i in range(N):
        K1X = Xas(Y[i],Z[i])
        K1Y = Yas(X[i],Y[i],a)
        K1Z = Zas(X[i],Z[i],b,r)

        K2X = Xas(Y[i]+(h/2)*K1Y,Z[i]+(h/2)*K1Z)
        K2Y = Yas(X[i]+(h/2)*K1X,Y[i]+(h/2)*K1Y,a)
        K2Z = Zas(X[i]+(h/2)*K1X,Z[i]+(h/2)*K1Z,b,r)

        K3X = Xas(Y[i]+(h/2)*K2Y,Z[i]+(h/2)*K2Z)
        K3Y = Yas(X[i]+(h/2)*K2X,Y[i]+(h/2)*K2Y,a)
        K3Z = Zas(X[i]+(h/2)*K2X,Z[i]+(h/2)*K2Z,b,r)

        K4X = Xas(Y[i]+(h/2)*K3Y,Z[i]+(h/2)*K3Z)
        K4Y = Yas(X[i]+(h/2)*K3X,Y[i]+(h/2)*K3Y,a)
        K4Z = Zas(X[i]+(h/2)*K3X,Z[i]+(h/2)*K3Z,b,r)

        X[i+1] = X[i]+(h/6)*(K1X+(2*K2X)+(2*K3X)+K4X)
        Y[i+1] = Y[i]+(h/6)*(K1Y+(2*K2Y)+(2*K3Y)+K4Y)
        Z[i+1] = Z[i]+(h/6)*(K1Z+(2*K2Z)+(2*K3Z)+K4Z)
    
    return X,Y,Z

#Решение Генераторв
def eiler_GKPR(Xx,Yy,Zz,N,h,u,g):
    X,Y,Z = np.empty(N+1), np.empty(N+1),np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()
    Z[0] = np.random.rand()
    
    for i in range(N):
        X[i+1] = X[i] + h * Xx(X[i],Y[i],Z[i],u,g)
        Y[i+1] = Y[i] + h * Yy(X[i])
        Z[i+1] = Z[i] + h * Zz(X[i],Z[i])



    return X,Y,Z

def R_K_GKPR(Xas,Yas,Zas,N,h,u,g):
    X,Y,Z= np.empty(N+1), np.empty(N+1),np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()
    Z[0] = np.random.rand()

    for i in range(N):
        K1X = Xas(X[i],Y[i],Z[i],u,g)
        K1Y = Yas(X[i])
        K1Z = Zas(X[i],Z[i] )

        K2X = Xas(X[i] + (h/2)*K1X,Y[i] + (h/2) * K1Y,Z[i] + (h/2) * K1Z , u, g )
        K2Y = Yas(X[i]+(h/2)*K1X)
        K2Z = Zas(X[i]+(h/2)*K1X,Z[i]+(h/2)*K1Z)

        K3X = Xas(X[i] + (h/2)*K2X,Y[i] + (h/2) * K2Y,Z[i] + (h/2) * K2Z , u, g )
        K3Y = Yas(X[i]+(h/2)*K2X)
        K3Z = Zas(X[i]+(h/2)*K2X,Z[i]+(h/2)*K2Z)

        K4X = Xas(X[i] + (h/2)*K3X,Y[i] + (h/2) * K3Y,Z[i] + (h/2) * K3Z , u, g )
        K4Y = Yas(X[i]+(h/2)*K3X)
        K4Z = Zas(X[i]+(h/2)*K3X,Z[i]+(h/2)*K3Z)

        X[i+1] = X[i]+(h/6)*(K1X+(2*K2X)+(2*K3X)+K4X)
        Y[i+1] = Y[i]+(h/6)*(K1Y+(2*K2Y)+(2*K3Y)+K4Y)
        Z[i+1] =  Z[i]+(h/6)*(K1Z+(2*K2Z)+(2*K3Z)+K4Z)
    return X,Y,Z

#Обратный метод Эйлера для систем

def O_eiler_for_Ueda(F,G,E,N,A,w,h):
    X,Y = np.empty(N+1), np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()
    for i in range(N):


        nX = X[i]
        NX = 2*nX
        ny = Y[i]
        Ny = 2*ny


        while(  (abs(nX -NX) > E) and (abs(ny - Ny) > E) ):
            NX = nX
            Ny = ny

            nX = X[i]+h * F(Ny)
            ny = Y[i] + h *G(NX,Ny,A,w,h*i)

        X[i+1] = nX
        Y[i+1] = ny
    
    return X,Y

def O_eiler_for_ressler(E,N,h,a,r,b=0.2):
    X,Y,Z = np.empty(N+1), np.empty(N+1),np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()
    Z[0] = np.random.rand()
    for i in range(N):
        nX = X[i]
        NX = 2*nX

        ny = Y[i]
        Ny = 2*ny

        nz = Z[i]
        Nz = 2*nz

        while(  (abs(nX -NX) > E) and (abs(ny - Ny) > E) and ( abs(nz - Nz) > E ) ):
            NX = nX
            Ny = ny
            Nz = nz

            nX =  X[i]+ h * ( -Y[i] - Z[i] )

            ny = Y[i]+ h* ( X[i] + a * Y[i] )

            nz = Z[i] + h * (  b + (X[i] - r)*Z[i] )


        X[i+1] = nX
        Y[i+1] = ny
        Z[i+1] = nz
    return X,Y,Z

def O_eiler_for_GKPR(Xx,Yy,Zz,E,N,h,u,g):
    X,Y,Z = np.empty(N+1), np.empty(N+1),np.empty(N+1)
    X[0] = np.random.rand()
    Y[0] = np.random.rand()
    Z[0] = np.random.rand()
    for i in range(N):
        nX = X[i]
        NX = 2*nX
        ny = Y[i]
        Ny = 2*ny
        nz = Z[i]
        Nz = 2*nz
        while(  (abs(nX -NX) > E) and (abs(ny - Ny) > E) and ( abs(nz - Nz) > E ) ):
            NX = nX
            Ny = ny
            Nz = nz
            nX = X[i] + h * Xx(X[i],Y[i],Z[i],u,g)
            ny = Y[i] + h * Yy(X[i])
            nz = Z[i] + h * Zz(X[i],Z[i])
        X[i+1] = nX
        Y[i+1] = ny
        Z[i+1] = nz
    return X,Y,Z


    X = np.empty(N)
    X[0] = 0.1

    for i in range(1,N):

        nX = X[i-1]
        NX = 2*nX

        while(  (abs(nX -NX) >= E) ):
            NX = nX
            nX = X[i-1]+h * F(nX)
  
        X[i] = nX
        NX = 0

    return X

    X= np.empty(N+1)
    X[0] = 0.1

    for i in range(0,N):

        X[i+1] = X[i]+h*F(X[i])

    return X






#Общие параметры
h  = 0.001 #шаг
E  = 0.000001 # погрешность для обратного метода
N = 100000 #число точек

#Решение систем Рёсслера
#параметры
a = 0.2
r = 4
#D,B,K =  eiler_Ressler(N,h,a,r)
#D,B,K = Runge_kutte_for_ressler(Xaxs,Yaxs,Zaxs,N,h,a,r)
#D,B,K = O_eiler_for_ressler(E, N, h,a,r)

#Решение генератора
#Параметры
u = 0.15
g = 0.93
#D,B,K = eiler_GKPR(Xax,Yax,Zax,N,h,u,g)
#D,B,K = R_K_GKPR(Xax,Yax,Zax,N,h,u,g)
#D,B,K = O_eiler_for_GKPR(Xax,Yax,Zax,E,N,h,u,g)

#График для трёхмерных систем
#Fig = plt.figure()
#Ax = Fig.add_subplot(111, projection='3d')
#Ax.plot(D[1000:],B[1000:],K[1000:])

#Решение Уеды
#Параметры
a = 420
w = 2.2
#D,B = eiler_ueda(F_eiler,G_Eiler,N,a,w,h)
#D,B = Runge_kutte_for_Ueda(F_eiler,G_Eiler,N,a,w,h)
#D,B = O_eiler_for_Ueda(F_eiler,G_Eiler,E,N,a,w,h)

#plt.plot(D[1000:],B[1000:])
plt.show()


