import pylab as pyl
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

c=3e8
e=1.6*1e-19
b=1

# bellurie
fig, ax1= plt.subplots()
ax1.set_title(r"Campo elettrico di carica in moto")
ax1.set_xlabel('t [s] ')
ax1.set_ylabel('E [V/m]')

# quantità preliminari

def gamma(beta):
    return 1/np.sqrt(1-beta**2)
def R(t,beta, b):
    return (gamma(beta)**2*(beta*c)**2*t**2 + b**2)**(3/2)

# Campi
def E_x(t,beta,b):
    return -e*gamma(beta)*beta*c*t/R(t,beta,b)
def E_y(t,beta,b):
    return e*b*gamma(beta)/R(t,beta,b)

# vettori utili
tt=np.linspace(-1e-7,1e-7,2000)
beta=np.linspace(0.1,0.9, 3)
print(beta)

pyl.plot(tt,E_x(tt,beta[0],b), color='blue',  linewidth = '2.5' , label=r"$E_x$ con $\beta = 0.1$")
pyl.plot(tt,E_y(tt,beta[0],b), color='black',  linewidth = '2.5' , label=r"$E_y$ con $\beta = 0.1$")

pyl.plot(tt,E_x(tt,beta[1],b), color='blue',  linewidth = '1' , label=r"$E_x$ con $\beta = 0.5$")
pyl.plot(tt,E_y(tt,beta[1],b), color='black',  linewidth = '1' , label=r"$E_y$ con $\beta = 0.5$")

pyl.plot(tt,E_x(tt,beta[2],b), color='blue',  linewidth = '0.5' , label=r"$E_x$ con $\beta = 0.9$")
pyl.plot(tt,E_y(tt,beta[2],b), color='black',  linewidth = '0.5' , label=r"$E_y$ con $\beta = 0.9$")

pyl.legend(loc = 'upper right')

plt.show()
