import difmetods as dif
import numpy as np
import matplotlib.pyplot as plt
################################################################
#генерирование данных
sins = [np.sin(i) for i in np.arange(-5,5,0.01)]
cos = [np.cos(i) for i in np.arange(-5,5,0.01)]
msins = [-np.sin(i) for i in np.arange(-5,5,0.01)]
T1 = [i for i in np.arange(-5,5,0.01)]
#################################################################
#Константы
m1 = 5
m21 = int((m1-1)/2)
AorB = 1
##################################################################

dsin = dif.parabdif2(sins,T1,m1,2)
d2sin = dif.parabdif2(sins,T1,m1,AorB)

plt.plot(T1[m21:len(T1)-m21], dsin,color = 'red')
plt.plot(T1,cos,"--",color = 'green')
plt.plot(T1,sins)
plt.plot(T1[m21:len(T1)-m21],d2sin,color = 'blue')
plt.plot(T1,msins,'--')
plt.show()