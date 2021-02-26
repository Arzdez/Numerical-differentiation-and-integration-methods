import math
import numpy as np
import metod

def f(x): return math.sin(x)
def df(x): return math.cos(x)



T = np.arange(-6, 6, 0.01)
K = np.arange(0,1,0.01)
print(metod.sravnen(f, df, T, K, M=11))




