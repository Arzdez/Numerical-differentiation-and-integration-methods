import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams
import numpy as np


def linedif(f, df, T, m=5, Grid=False, Gause=0, Eq=True):

    def M2(m): return int((m - 1) // 2)

    def Kn(X, n, z):
        M = [i for i in range(-m2, m2 + 1)]
        return (1 / z) * sum([X[int(i + n)] * i for i in M])

    def KnNoEq(n, T, X):

        m2 = M2(m)
        S1, S2, S3, S4 = 0, 0, 0, 0

        for j in range(n - m2, n + m2 + 1):
            S1 += X[j] * T[j]
            S2 += X[j]
            S3 += T[j]
            S4 += T[j] ** 2

        return (S1 - (1 / m) * S2 * S3) / (S4 - ((1 / m) * (S3 ** 2)))

    # Константы
    L = len(T)
    dt = T[1] - T[0]
    m2 = M2(m)
    Z = (2 * dt) * sum(i ** 2 for i in [Mi for Mi in range(1, m2 + 1)])

    # Получаем ряд значений и их сумма с шумом
    X = [f(x) for x in T]
    X_Noise = [i + np.random.normal(0, Gause, 1) for i in X]

    # Получаем производную
    dT = [T[i] for i in range(m2, L - m2)]

    # Получаем диф. ряд
    if (Eq): dX = [Kn(X_Noise, n, Z) for n in range(m2, len(dT) + m2)]
    else: dX = [KnNoEq(n, T, X_Noise) for n in range(m2, len(dT) + m2)]

    # Идеальная производная
    DX = [df(x) for x in T]

    # Считаем квадратичную ошибку
    Kfn = sum([(DX[i + m2] - dX[i]) ** 2 for i in range(len(dX))])

    # Рисуем
    TITLE = "-"
    if (Grid): PloterGraphFdF(TITLE,T,X,dT,dX,T,DX,Kfn)
    return Kfn


def LeftDiff(f, df, T, Grid=False, Gause=0):

    # Константы
    L = len(T)

    # Получаем ряд значений и их сумма с шумом
    X = [f(x) for x in T]
    X_Noise = [i + np.random.normal(0, Gause, 1) for i in X]

    # Получаем производную
    dT = [T[i] for i in range(1, L)]
    dX = [(X_Noise[i + 1] - X_Noise[i]) / (T[i + 1] - T[i]) for i in range(len(dT))]

    # Идеальная производная
    DX = [df(x) for x in T]

    # Считаем квадратичную ошибку
    Kfn = sum([(DX[i + 1] - dX[i]) ** 2 for i in range(len(dX))])

    # Рисуем
    if (Grid): PloterGraphFdF("LeftDiv",T,X,dT,dX,T,DX,Kfn)
    return Kfn


def RightDiff(f, df, T, Grid=False, Gause=0):

    # Константы
    L = len(T)

    # Получаем ряд значений и их сумма с шумом
    X = [f(x) for x in T]
    X_Noise = [i + np.random.normal(0, Gause, 1) for i in X]

    # Получаем производную
    dT = [T[i] for i in range(1, L - 1)]
    dX = [(X_Noise[i] - X_Noise[i - 1]) / (T[i] - T[i - 1]) for i in range(1, L - 1)]

    # Идеальная производная
    DX = [df(x) for x in T]

    # Считаем квадратичную ошибку
    Kfn = sum([(DX[i + 1] - dX[i]) ** 2 for i in range(len(dX))])

    # Рисуем
    if (Grid): PloterGraphFdF("RightDiv",T,X,dT,dX,T,DX,Kfn)
    return Kfn


def CentralDiff(f, df, T, Grid=False, Gause=0):

    # Константы
    L = len(T)

    # Получаем ряд значений и их сумма с шумом
    X = [f(x) for x in T]
    X_Noise = [i + np.random.normal(0, Gause, 1) for i in X]

    # Получаем производную
    dT = [T[i] for i in range(1, L - 1)]
    dX = [(X_Noise[i + 1] - X_Noise[i - 1]) / (T[i + 1] - T[i - 1]) for i in range(1, L - 1)]

    # Идеальная производная
    DX = [df(x) for x in T]

    # Считаем квадратичную ошибку
    Kfn = sum([(DX[i + 1] - dX[i]) ** 2 for i in range(len(dT))])

    # Рисуем
    if (Grid): PloterGraphFdF("CentralDiv",T,X,dT,dX,T,DX,Kfn)
    return Kfn


def sravnen(f, df, T, K, M=3, Eq=False):
    LD, RD, CD, DD = [], [], [], []
    for t in K:
        LD.append(LeftDiff(f, df, T=T, Gause=t))
        #RD.append(RightDiff(f, df, T=T, Gause=t))
        CD.append(CentralDiff(f, df, T=T, Gause=t))
        DD.append(linedif(f, df, T=T, m=M, Gause=t, Eq=Eq))

    # Рисуем
    plt.title("Сравнение точности методов в зависимости от уровня шума")
    plt.ylabel("df - f'")
    plt.xlabel("Коэф. шума")
    plt.plot(K, LD, '--', label="Слева",color = 'green')
    #plt.plot(K, RD, '--', label="Right Div")
    plt.plot(K, CD, '--', label="симметричная",color = 'blue')
    plt.plot(K, DD, '--', label="с помощью линейной функции",color = 'red')
    plt.legend()
    plt.show()


def PloterGraphFdF(Title, X1, Y1, X2, Y2, X3, Y3, Kfn):
    plt.title(Title)
    plt.xlabel("Kfn: " + str(Kfn))
    plt.plot(X1, Y1, color="green", label="f(x)")
    plt.plot(X2, Y2, color="red", label="df(x)")
    plt.plot(X3, Y3, "--", label="f'(x)")
    plt.legend()
    plt.show()


# def DarkMode():
#     rcParams['figure.edgecolor'] = "333333"
#     rcParams['figure.facecolor'] = "333333"
#     rcParams['figure.figsize'] = 15, 9

#     rcParams['text.color'] = "CCCC00"

#     rcParams['axes.labelcolor'] = "CCCC00"
#     rcParams['axes.edgecolor'] = "ffffff"
#     rcParams['axes.facecolor'] = "222222"

#     rcParams['savefig.edgecolor'] = "222222"
#     rcParams['savefig.facecolor'] = "222222"

#     rcParams['xtick.color'] = "CCAA00"
#     rcParams['ytick.color'] = "CCAA00"

#     rcParams['xtick.minor.visible'] = True
#     rcParams['ytick.minor.visible'] = True

#     rcParams['boxplot.meanline'] = True
#     rcParams['figure.frameon'] = False
#     rcParams['grid.color'] = "055212"

#     rcParams['font.size'] = 16
