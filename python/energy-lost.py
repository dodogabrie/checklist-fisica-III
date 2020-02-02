import pylab as pyl
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#mev
me = 0.510 #9.1e-31 kg
mu = 931.5 #1.66e-27 kg
hbar = 1.05*1e-34 # J*s
c = 3e8 # m/s
e = 1.6*1e-19 # C
Na = 6.02*1e23
alpha = 1/137
re = 2.82*1e-15

# caso di particelle con energia inferiore a soglia:
def NoSoglia(A, Z, L, z, a, E):
    gamma = E/(me)
    n = Na/A
    cost =n*Z**2*z**2/(a**2)*alpha*re**2
    return 16/3*cost*np.log(gamma)*gamma*A*mu*L

# Piombo
Apb = 207
Zpb = 82
rhopb = 11.35*1e3 #kg/m**3
L = 2e-3 
Pb = np.array([Apb, Zpb, L])

# e- da 3.5 MeV
E = 3.5
print(NoSoglia(*Pb,1,1, E))
















