import pylab as pyl
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#proiettile
Ap = 1
Zp = 1
Dp = 0
#bersaglio
Ab = 16
Zb = 8
Db = -4.74

mp = 938.2
mu = 931.49
e0 = 8.85*1e-12
e = 1.6*1e-19
Q = 16.8

def r(A):
    if A == 1:
        return 1.25*1e-15
    else:
        return (1.25*A**(1/3) + 2)*1e-15  
def M(A,D):
    return A*mu + D
def E(z,Z,Ap,Ab):
    return (1 + M(Ap,Dp)/M(Ab,Db))*z*Z*e/(4*np.pi*e0*(r(Ap)+r(Ab)))*1e-6 + Q
print(E(Zp,Zb,Ap,Ab))
