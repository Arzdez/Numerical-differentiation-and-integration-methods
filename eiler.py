import math
import numpy as np
import matplotlib.pyplot as plt
import random
from math import log

from numpy.core.shape_base import atleast_1d
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
#Осцилятор Рёсслера 2
def Xaxs1(y,z,Xs,X1s,k):
    return (- y - z)- k*(X1s-Xs)

def Yaxs1(x,y,a):
    return ( x +a*y)

def Zaxs1(x,z,b,r):
    return ( b+(x-r)*z)

#Генератор Кияшко-Пиковского-Рабиновича
def Ff(z):
    return ( (8.592*z)-(22*z*z) + (14.408*(z*z*z)) )

def Xax(x,y,z,u,g):
    return (2*u*x + y - g*z) 

def Yax(x):
    return -x

def Zax(x,z):
    b = ((x - Ff(z) ) / 0.2 ) 
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


#Решение связанных систем
def eiler_sys(N,h, a, a1,r,r1,k,t,b=0.2,b1=0.2):

    X,Y,Z,X1,Y1,Z1 = np.empty(N+1), np.empty(N+1),np.empty(N+1),np.empty(N+1), np.empty(N+1),np.empty(N+1)
    for i in range(t):
        X[i] = np.random.rand()
        Y[i] = np.random.rand()
        Z[i] = np.random.rand()
        X1[i] = np.random.rand()
        Y1[i] = np.random.rand()
        Z1[i] = np.random.rand()

    for i in range(t,N):
        X[i+1] = X[i]+ h * ( -Y[i] - Z[i] )
        Y[i+1] = Y[i]+ h* ( X[i] + a * Y[i] )
        Z[i+1] = Z[i] + h * (  b + (X[i] - r)*Z[i] )

        X1[i+1] = X1[i] + h * ( ( -Y1[i] - Z1[i] ) - (X[i-t]*k ) )
        Y1[i+1] = Y1[i]+ h * ( X1[i] + a1 * Y1[i] )
        Z1[i+1] = Z1[i] + h * (  b1 + (X1[i] - r1)*Z1[i] )
    
    return X,Y,Z,X1,Y1,Z1

def RK_sys(Xas,Yas,Zas,Xas1,Yas1,Zas1,N,h,a,a1,r,r1,k,t,b=0.2,b1=0.2):

    X,Y,Z,X1,Y1,Z1= np.empty(N+1), np.empty(N+1),np.empty(N+1) ,np.empty(N+1), np.empty(N+1),np.empty(N+1)
    for i in range(t):
        X[i] = np.random.rand()
        Y[i] = np.random.rand()
        Z[i] = np.random.rand()
        X1[i] =np.random.rand()
        Y1[i] =np.random.rand()
        Z1[i] =np.random.rand()

    for i in range(t,N):
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

        K1X1 = Xas1( Y1[i], Z1[i], X[i-t], X1[i] , k)
        K1Y1 = Yas1( X1[i] , Y1[i], a1)
        K1Z1 = Zas1( X1[i] , Z1[i], b1 , r1)

        K2X1 = Xas1( Y1[i]+(h/2) * K1Y1 , Z1[i] + (h/2)*K1Z1 , X[i-t]+(h/2) * K1X, X1[i]+(h/2) * K1X1 , k)
        K2Y1 = Yas1( X1[i]+(h/2) * K1X1,Y1[i] + (h/2) * K1Y1, a1)
        K2Z1 = Zas1( X1[i]+(h/2) * K1X1,Z1[i] + (h/2) * K1Z1, b1 , r1)

        K3X1 = Xas1( Y1[i]+(h/2)  * K2Y1,Z1[i] + (h/2)*K2Z1 , X[i-t]+(h/2) * K2X,X1[i]+(h/2) * K2X1 , k)
        K3Y1 = Yas1(X1[i]+(h/2)*K2X1,Y1[i]+(h/2)*K2Y1,a1)
        K3Z1 = Zas1(X1[i]+(h/2)*K2X1,Z1[i]+(h/2)*K2Z1,b1,r1) 

        K4X1 = Xas1( Y1[i]+(h/2)*K3Y1 , Z1[i]+(h/2)*K3Z1 , X[i-t]+(h/2)*K3X ,X1[i]+(h/2) * K3X1 , k )
        K4Y1 = Yas1( X1[i]+(h/2)*K3X1 , Y1[i]+(h/2)*K3Y1 , a1 )
        K4Z1 = Zas1( X1[i]+(h/2)*K3X1 , Z1[i]+(h/2)*K3Z1, b1, r1)
        
        X1[i+1] = X1[i]+(h/6)*(K1X1+(2*K2X1)+(2*K3X1)+K4X1)
        Y1[i+1] = Y1[i]+(h/6)*(K1Y1+(2*K2Y1)+(2*K3Y1)+K4Y1)
        Z1[i+1] = Z1[i]+(h/6)*(K1Z1+(2*K2Z1)+(2*K3Z1)+K4Z1)
        
        

    return X,Y,Z, X1,Y1,Z1
#Обратный метод Эйлера для систем

#def O_eiler_for_Ueda(F,G,E,N,A,w,h):
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
#
#def O_eiler_for_ressler(E,N,h,a,r,b=0.2):
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
#
#def O_eiler_for_GKPR(Xx,Yy,Zz,E,N,h,u,g):
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
#





#Общие параметры











#D,B,K = Runge_kutte_for_ressler(Xaxs,Yaxs,Zaxs,N,h,a,r)
#D,B,K = O_eiler_for_ressler(E, N, h,a,r)

#Решение генератора
#Параметры
#u = 0.14
#g = 0.9
#D,B,K = eiler_GKPR(Xax,Yax,Zax,N,h,u,g)
#D1,B1,K1 = R_K_GKPR(Xax,Yax,Zax,N,h,u,g)
#D,B,K = O_eiler_for_GKPR(Xax,Yax,Zax,E,N,h,u,g)

h  = 0.001 #шаг
N = 100000 #число точек

#1 система
a =  0.32
r =  4.9
#2 система
a1 = 0.2
r1 = 2
#D,B,K =  eiler_Ressler(N,h,a,r)
#D,B,K,D1,B1,K1 = eiler_sys(N,h,a,a1,r,r1,k)
D1,B1,K1,D11,B11,K11 = eiler_sys(N,h,a,a1,r,r1,0.2,0)
D12,B12,K12,D112,B112,K112 = eiler_sys(N,h,a,a1,r,r1,0.2,200)
D13,B13,K13,D113,B113,K113 = eiler_sys(N,h,a,a1,r,r1,0.2,300)
D14,B14,K14,D114,B114,K114 = eiler_sys(N,h,a,a1,r,r1,0.2,400)
D15,B15,K15,D115,B115,K115 = eiler_sys(N,h,a,a1,r,r1,0.2,900)

#D1,B1,K1,D11,B11,K11 = RK_sys(Xaxs,Yaxs,Zaxs,Xaxs1,Yaxs1,Zaxs1,N,h,a,a1,r,r1,1,0)
#
#D12,B12,K12,D112,B112,K112 = RK_sys(Xaxs,Yaxs,Zaxs,Xaxs1,Yaxs1,Zaxs1,N,h,a,a1,r,r1,1,100)
#
#D13,B13,K13,D113,B113,K113 = RK_sys(Xaxs,Yaxs,Zaxs,Xaxs1,Yaxs1,Zaxs1,N,h,a,a1,r,r1,1,200)
#
#D14,B14,K14,D114,B114,K114 = RK_sys(Xaxs,Yaxs,Zaxs,Xaxs1,Yaxs1,Zaxs1,N,h,a,a1,r,r1,1,300)
#
#D15,B15,K15,D115,B115,K115 = RK_sys(Xaxs,Yaxs,Zaxs,Xaxs1,Yaxs1,Zaxs1,N,h,a,a1,r,r1,1,1000)
#График для трёхмерных систем
Fig = plt.figure()
Ax = Fig.add_subplot(251, projection='3d')
Az = Fig.add_subplot(256, projection='3d')
Az1 = Fig.add_subplot(257, projection='3d')
Az2 = Fig.add_subplot(258, projection='3d')
Az3 = Fig.add_subplot(259, projection='3d')
Az4 = Fig.add_subplot(2,5,10, projection='3d')

Ax.plot(D1[30000:],B1[30000:],K1[30000:])
Az.plot(D11[30000:],B11[30000:],K11[30000:])
Az1.plot(D112[30000:],B112[30000:],K112[30000:])
Az2.plot(D113[30000:],B113[30000:],K113[30000:])
Az3.plot(D114[30000:],B114[30000:],K114[30000:])
Az4.plot(D115[30000:],B115[30000:],K115[30000:])
    
#Решение Уеды
#Параметры
a = 420
w = 2.2
#D,B = eiler_ueda(F_eiler,G_Eiler,N,a,w,h)
#D,B = Runge_kutte_for_Ueda(F_eiler,G_Eiler,N,a,w,h)
#D,B = O_eiler_for_Ueda(F_eiler,G_Eiler,E,N,a,w,h)

#plt.plot(D[1000:],B[1000:])
plt.show()


