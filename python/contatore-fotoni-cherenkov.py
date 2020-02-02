import pylab as pyl
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

beta = 0.99
n=1.5
l1 = 300e-9
l2 = 600e-9
L = 1e-2
P = 0.3
def Count(L,l1,l2,P,beta, n):
    return 2*np.pi*1/137*( 1/l1 - 1/l2 )*(1- 1/(beta**2*n**2))*P

print(Count(L,l1,l2,P,beta, n))
