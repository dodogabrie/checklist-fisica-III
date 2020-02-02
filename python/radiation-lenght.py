import pylab as pyl
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

Ap = 207
Zp = 82
As = 28
Zs = 14


NA = 6.02*1e23
re = 2.82*1e-13
alpha = 1/137

def rhoX0(A,Z):
    return 1/(16/3*NA/(A)*(Z)**2*alpha*(re)**2*np.log(192/(Z**(1/3))))
print(rhoX0(As,Zs))
