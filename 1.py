import difmetods as dif
import numpy as np
import matplotlib.pyplot as plt
sins = []
cos = []
msins = []
T1 = []
m1 = 5
m21 = int((m1-1)/2)
for i in np.arange(-5,5,0.01):
    sins.append(np.sin(i))
    cos.append(np.cos(i))
    msins.append(-np.sin(i))
    T1.append(i)
dsin = dif.parabdif2(sins,T1,m1,2)
d2sin = dif.parabdif2(dsin,T1[m21:len(T1)-m21],m1,2)

plt.plot(T1[m21:len(T1)-m21], dsin,color = 'red')
plt.plot(T1,cos,"--",color = 'green')
plt.plot(T1,sins)
plt.plot(T1[2*m21:len(T1)-2*m21],d2sin,color = 'blue')
plt.plot(T1,msins,'--')
plt.show()