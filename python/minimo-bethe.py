import pylab as pyl
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


Z1 = 82
A1 = 207

Z2 = 14
A2 = 28

Z3 = 7
A3 = 14

def beth(A,Z):
    return 0.307*Z/A*(np.log(7.8*1e5* 1/(Z**(0.9))) - 1)
print(beth(A2,Z2))
