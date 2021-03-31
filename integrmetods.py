import numpy as np
import matplotlib.pyplot as plt
import random
#служебные функции для метода Гауса( содержат коэфициенты )
def xi(n):
    if n == 1:
        return [0]
    if n == 2:
        return [-0.5773503,0.55773503]
    if n == 3:
        return [-0.7745967,0,0.7745967]
    if n == 4:
        return [-0.8611363,-0.3399810,0.3399810,0.8611363]
    if n == 5:
        return [-0.9061798,-0.5384693,0,0.5384693,0.9061798]
def cin2(n):
    if n == 1:
        return [2]
    if n == 2:
        return [1,1]
    if n == 3:
        return [0.5555556,0.8888889,0.5555556]
    if n == 4:
        return [0.3478548,0.6521451,0.6521451,0.3478548]
    if n == 5:
        return [0.4786287,0.2369269,0.5688888,0.2369269,0.4786287]
#служебная функия для метода N-K(содержит коэфициенты)
def cin(n):
    if n == 0:
        return [1]
    if n == 1:
        return [1,1]
    if n == 2:
        return [1,4,1]
    if n == 3:
        return [1,3,3,1]
    if n == 4:
        return [7,32,12,32,7]
    if n == 5:
        return [19,75,50,50,75,19]
#служембная функция для метода монте-Карло
def get_rand_number(min_value, max_value):
    range = max_value - min_value
    choice = random.uniform(0,1)
    return min_value + range*choice
# правые прямоугольники
def Iright(Fx,a,b,n):
    dx=(b-a)/n
    result = 0
    for i in range(n):
        result += Fx(dx*i)

    return(result*dx)
#Центральные прямоугольники
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
#метод трапеций
def Itrapeciya(Fx,a,b,n):
    dx = (b-a)/n
    sumy = 0
    for i in range(n):
        sumy+=Fx(dx*i)
    I = dx*( ((Fx(a)+Fx(b))/2)+sumy)
    return( I)
#Метод Симпсона
def Isimpson(Fx,a,b,n): 
    h=((b-a)/n)
    s=Fx(a)+Fx(b)
    for i in range(1,n):
        if i%2==0:
            s+= 2*Fx(a+h*i) 
        else:
            s+= 4*Fx(a+h*i)
    return(s*h/3)
#Монте Карло первый способом
def monte_carlo(Fx,a,b,num_samples):
    sum_of_samples = 0
    for i in range(num_samples):
        x = get_rand_number(a,b)
        sum_of_samples += Fx(x)
    
    return (b - a) * float(sum_of_samples/num_samples)
#Монте Карло прямоугольником
def monte_Karlo(Fx,a,b,num_sample):
    i1 = 0
    for i in range(num_sample):
        x = get_rand_number(a,b)
        y = get_rand_number(0,b)
        if y<=Fx(x):
            i1+=1       
    return(b-a)*b*(i1/num_sample)
#Ньютон котес
def Nuton_Kotes(Fx,a,b,N,dig):
    suma = 0
    C = sum(cin(dig))
    h = (b - a) / (N*dig)
    mlp = dig/C * h

    for j in range(0,N):
        partsum = 0
        for i in range(0,dig+1):
            partsum += cin(dig)[i] * Fx(a + ( i + j* dig ) * h)
        suma+=partsum 

    return suma * mlp
#Гаус
def Gaus(Fx,a,b,N,dig):
    suma = 0
    mlp = (b-a)/(2*N)
    for j in range(0,N):
        partsum = 0
        for i in range(0,dig):
            xd = a + ( ((xi(dig)[i]+1)*(b-a))/2 )
            partsum += cin2(dig)[i] * Fx( xd)
        
        suma+=partsum 

    return( suma * mlp)



if __name__ == "__main__":

    def sins(x):
        return np.sin(x)
    a = 0
    b = np.pi/2
    graph = []
    x = []
    for i in range(100):
        graph.append(sins(i))
        x.append(i)

    print(Iright(sins,a,b,1000))
    print(Imidle(sins,a,b,1000))
    print(Itrapeciya(sins,a,b,1000))
    print(Isimpson(sins,a,b,1000))
    print(monte_carlo(sins,a,b,1000))
    print(monte_Karlo(sins,a,b,1000))
    print(Nuton_Kotes(sins,a,b,1000,3))
    print(Gaus(sins,a,b,1000,3))