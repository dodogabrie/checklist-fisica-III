import pylab as pyl
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# velocità della luce
c = 3*1e8
# costante di plank ridotta
h = 6.6*1e-16 # eV*s 
# costante dielettrica nel vuoto
e0 = 8.85*1e-12
# carica elementare
e = 1.6*1e-19
# elettrone
z = 1
T = 1e3 
m = 0.511 # MeV/c**2
# Piombo
Z = 82
A = 207

def beta(T, m):
    E = T + m
    return np.sqrt(1-(m/E)**2)
def sigma_th(z,Z,T, theta):
    return (z*Z*e/(4*np.pi*e0))**2*(1/(4*T*1e6))**2*1/(np.sin(theta/2))**4
def sigma_mott(z,Z,T,theta, m):
    return sigma_th(z,Z,T,theta)*(1-beta(T,m)**2*np.sin(theta/2)**2)
# funzione del raggio nucleare
def R(A):
    return (1.25*A**(1/3)+2)*1e-15
# funzione del fattore di forma
def F(theta, A, m, T):
    p = np.sqrt((T+m)**2 - m**2)*1e6 #eV/c
    x = 2*R(A)*p/(h*c)*np.sin(theta/2) # adim
    return 3*(np.sin(x)/x**3-np.cos(x)/x**2)
def sigma_tot(theta, z, Z, A, m, T):
    return sigma_mott(z,Z,T,theta,m)*np.abs(Z*F(theta, A, m, T))**2

# bellurie
pyl.title(r'Misura del raggio con fattore di forma')
plt.rc('font',size=10)
pyl.yscale('log')    #SCALA LOGARITMICA SOLO SU X
plt.ylim(1e-4, 1e4)
plt.xlim(0.4, 2)
plt.ylabel(r'$\frac{d\sigma}{d \Omega}$ [barn]')
plt.xlabel(r'$\theta$ [$^o$]')
pyl.minorticks_on()

theta = np.linspace(0.1/(2*np.pi), 3/(2*np.pi), 500000)


                                # per decommentare togli cancelletti rendendo pari qui-->|
#----- per vedere che il minimo è a 4.5--------------------------------------------------|
# bellurie                                                                              #|
#pyl.title(r'Fattore di forma per sfera uniformemente carica')                           #|
#plt.rc('font',size=10)                                                                  #|
#plt.ylim(-0.13,1.1)                                                                     #|
#plt.xlim(0.1, 14)                                                                       #|
#plt.ylabel(r'F(qR)')                                                                    #|
#plt.xlabel(r'qR')                                                                       #|
#pyl.minorticks_on()                                                                     #|
# "x = Rq"                                                                              #|
#x = np.linspace(0.1,15, 5000)                                                           #|
#xx = np.linspace(-1, 20, 3000)                                                          #|
#def F1(x):                                                                              #|   
#    return 3*(np.sin(x)/x**3-np.cos(x)/x**2)                                            #|
#pyl.plot(x, F1(x), color='black')                                                       #|
#pyl.plot(xx, np.zeros(len(xx)), color='black', linestyle = '--', linewidth = '0.6')     #|
#pyl.show()                                                                              #|
#----------------------------------------------------------------------------------------|


#----------- Plot della sezione Mott -----------------------------------------|
#pyl.plot(theta*2*np.pi, sigma_mott(z, Z, T, theta, m)*1e28, color='blue')    #|
#pyl.show()                                                                   #|     
#-----------------------------------------------------------------------------|

#----------- Plot Del fattore di forma -----------------------------------|
#pyl.plot(theta*2*np.pi, F(theta,A, m, T), color='blue')                  #|
#pyl.plot(theta*2*np.pi, np.zeros(len(theta)), color='black')             #|
#pyl.show()                                                               #|   
#-------------------------------------------------------------------------|

#----------- Plot della sezione tot --------------------------------------------|
pyl.plot(theta*2*np.pi, sigma_tot(theta, z, Z, A, m, T)*1e28, color='blue')    #|
pyl.show()                                                                     #|
#-------------------------------------------------------------------------------|


a = 4.5
p = np.sqrt((T+m)**2 - m**2)*1e6 #eV/c
thh = np.arcsin(a*h*c/(2*R(A)*p))
print("angolo minimo:")
print(thh*2*np.pi*2)
