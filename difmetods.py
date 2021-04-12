import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tkinter.filedialog import askopenfilename
from scipy import interpolate
#i в Ti и В Xi равняется степени
#n - отрезок
#m - окно
#T - массив временных отрезков
#X - массив значений функции
#BorA = выбор режима работы функции просчёта производной с помощью параболы:
#  1 - оценка второй производной, 2 - оценка первой производной

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
        #расчёт коэффициента A
        A1 = XX[2] - om * TT[2] * XX[0]
        A2 = TT[3] - om * TT[2] * TT[1]
        A3 = XX[1] - om * XX[0] * TT[1]
        A4 = TT[2] - om * (TT[1] ** 2)
        A5 = TT[4] - om * (TT[2] ** 2)
        A6 = TT[3] - om * TT[2] * TT[1]
        A7 = om * TT[2] * TT[1] - TT[3]
        A8 = TT[2] - (om * (TT[1] ** 2))

        an = ((A1 - (A2 * (A3 / A4))) / (A5 + (A6 * (A7 / A8))))
        #Расчёт коэфициента B
        B1 = XX[1] - om * XX[0] * TT[1]
        B2 = TT[2] - om * (TT[1] ** 2 )
        B3 = om * TT[2] * TT[1] - TT[3]
        B4 = TT[2] - om * (TT[1] ** 2)
        bn = (B1 / B2) + (an * (B3 / B4))

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
    ################################################################
    #генерирование данных
    #sins = [np.sin(i) for i in np.arange(-5,5,0.01)]
    #cos = [np.cos(i) for i in np.arange(-5,5,0.01)]
    #msins = [-np.sin(i) for i in np.arange(-5,5,0.01)]
    #T1 = [i for i in np.arange(-5,5,0.01)]
    PATH = askopenfilename()

    FILES = pd.read_excel(PATH, sheet_name="Отчет")

    X = [float(FILES['Unnamed: 0'][i]) for i in range(10, len(FILES['Unnamed: 0']))]
    Y = [float(FILES['Unnamed: 1'][i]) for i in range(10, len(FILES['Unnamed: 1']))]
    #################################################################
    InX = np.linspace(min(X), max(X), 2000)                                 # Генерируем ряд, по которому произведём интерполяцию
    Tck = interpolate.splrep(X, Y)                                         # Заготовка для интерполяции
    InY = interpolate.splev(InX, Tck)                                      # Получаем интерполированные данные
    #Константы
    m1 = 5
    #m21 = int((m1-1)/2)
    AorB = 1
    ##################################################################

    dsin = parabdif2(InY,InX,m1,1)
    d2sin = parabdif2(InY,InX,m1,2)
    d22sin = parabdif2(d2sin,InX[2:len(InX)-2],m1,2)

    plt.plot(InX,InY)
    plt.plot(InX[2:len(InX)-2], d2sin,color = 'red')
    plt.plot(InX[4:len(InX)-4], d22sin)
    #plt.plot(T1,cos,"--",color = 'green')
    #plt.plot(T1,sins)
    #plt.plot(T1[m21:len(T1)-m21],d2sin,color = 'blue')
    #plt.plot(T1,msins,'--')
    plt.show()


