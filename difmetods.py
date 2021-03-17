import math
import numpy as np
import matplotlib.pyplot as plt
#i в Ti и В Xi равняется степени
#n - отрезок
#m - окно
#T - массив временных отрезков
#X - массив значений функции
#BorA = выбор режима работы функции просчёта 1 - оценка второй производной, 2 - оценка первой производной

##############################################################################
#Для аппроксимации параболой
def Ti(i,n,m2,T): return sum([T[j]**i for j in range(n - m2, n + m2+1)])
def Xi(i,n,m2,T,X): return sum([X[j]*(T[j]**i) for j in range(n - m2, n + m2+1)])
def BAn(n,om,m2,T,X,bora):
    TT = np.empty(5)
    XX = np.empty(5)

    for i in range(0,5):
        TT[i]=Ti(i,n,m2,T)
        XX[i]=Xi(i,n,m2,T,X)
    #просчёт А
    if bora == 1:
        K1 = XX[2] - om * TT[2] * XX[0]
        K2 = TT[3] - om * TT[2] * TT[1]
        K3 = XX[1] - om * XX[0] * TT[1]
        K4 = TT[2] - om * (TT[1] ** 2)
        K5 = TT[4] - om * (TT[2] ** 2)
        K6 = TT[3] - om * TT[2] * TT[1]
        K7 = (om * TT[2] * TT[1]) - TT[3]
        K8 = TT[2] - (om * (TT[1] ** 2))
        return ((K1 - (K2 * (K3 / K4))) / (K5 + (K6 * (K7 / K8))))*2
    #просчёт A и B
    else:
        K1 = XX[2] - om * TT[2] * XX[0]
        K2 = TT[3] - om * TT[2] * TT[1]
        K3 = XX[1] - om * XX[0] * TT[1]
        K4 = TT[2] - om * (TT[1] ** 2)
        K5 = TT[4] - om * (TT[2] ** 2)
        K6 = TT[3] - om * TT[2] * TT[1]
        K7 = om * TT[2] * TT[1] - TT[3]
        K8 = TT[2] - (om * (TT[1] ** 2))

        an = ((K1 - (K2 * (K3 / K4))) / (K5 + (K6 * (K7 / K8))))

        B1 = XX[1] - om*XX[0]*TT[1]
        B2 = TT[2] - om* (TT[1]**2)
        B3 = om*TT[2]*TT[1]-TT[3]
        B4 = TT[2] - om*(TT[1]**2)
        bn = (B1/B2) + (an*(B3/B4))

        return (2 * an * T[n] + bn)
def parabdif2(X,T,m,bora):
    m2 = int((m-1)/2)
    om = 1/m
    L = len(X)
    result = []
    for n in range(m2,L-m2):
        result.append(BAn(n,om,m2,T,X,bora))
    return result
##############################################################################
#Аппроксимация линейной функцией
def lindif(X,T,m):
    m2 = int((m-1)/2)
    k = []
    dt = T[1]-T[0]
    z = (2 * dt) * sum(i ** 2 for i in [Mi for Mi in range(1, m2 +1)])
    for n in range(m2,len(X)-m2):
        M = [i for i in range(-m2, m2+1 )]
        k.append ((1/z) * sum([X[i+n]*i for i in M]))
    return (k)  
###############################################################################
if __name__ == "__main__":
    print('huy')

    sins = []
    cos = []
    T1 = []
    m1 = 5
    m21 = int((m1-1)/2)
    msins = []
    for i in np.arange(-5,5,0.01):
        sins.append(np.sin(i))
        cos.append(np.cos(i))
        msins.append(-np.sin(i))
        T1.append(i)


