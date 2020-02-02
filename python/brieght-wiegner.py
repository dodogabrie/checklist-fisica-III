import pylab as pyl
import numpy as np
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

c = 3e8
e0 = 8.854e-12
e = 1.6e-19
me = 9.1e-31
re = 1/(4*np.pi*e0)* e**2/(me*c**2)
tau = (2/3)*re/c
L1 = 10e10
w0 = 1e15
L = tau*w0**2
barn = 1e30
lamda0 = 2*np.pi*c/w0

# thompson section
def sigma_th():
    return (8/3) *np.pi* re**2

# gamma_tot(w)
def Ltot(w):
    return L1 + w**2*L/(w0**2)

# aux. function anti-waste of time
def auxL(w):
    return 1/((w0**2 - w**2)**2 + w**2 * Ltot(w)**2)

# fit function
def sigma_tot(w):    
    return barn*4*np.pi*re*c*w**2*Ltot(w)*auxL(w)

def sigma_el(w):
    return barn*sigma_th()*w**4*auxL(w)

def sigma_abs(w):
    return barn*4*np.pi*re*w**2*auxL(w)*(c*Ltot(w) - 2/3*re*w**2 )


xx=np.logspace(14, 25 ,20000)
yy = np.logspace(0,15, 2000)
# bellurie
fig, ax1 = plt.subplots()
ax1.set_title(r"Bright-Wiegner")
#plt.rc('font',size=10)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim(1,4e14)
ax1.set_xlabel('pulse [Hz]')
ax1.set_xlim(1e14, 1.36e24)
x_label3 = (r"$\omega_0$", r"$\frac{1}{\tau}$")
x_tick = [w0, 1/tau]
plt.xticks(x_tick, x_label3)
pyl.plot(xx,sigma_tot(xx), color='red', linewidth = '1.4', label="$\sigma_{tot}$")
pyl.plot(xx,sigma_el(xx), color='blue',  linewidth = '1',label="$\sigma_{el}$")
pyl.plot(xx,sigma_abs(xx), color='black',  linewidth = '0.8',label="$\sigma_{abs}$")
pyl.legend(loc = 'upper right')
pyl.plot(xx, barn*np.ones(len(xx))*3/(2*np.pi)*lamda0**2*L/(L+L1), linestyle = '--', color='red',  linewidth = '0.5')
pyl.plot(xx, barn*np.ones(len(xx))*3/(2*np.pi)*lamda0**2*(L/(L+L1))**2, linestyle = '--', color='blue',  linewidth = '0.5')
pyl.plot(xx, barn*np.ones(len(xx))*sigma_th(), linestyle = '--', color='green',  linewidth = '0.5')
pyl.plot(np.ones(len(yy))*w0, yy , linestyle = '--', color='black',  linewidth = '0.6')
pyl.plot(np.ones(len(yy))*1/tau, yy , linestyle = '--', color='black',  linewidth = '0.6')


ax2 = ax1.secondary_yaxis('right')
#ax2 = plt.twinx()
ax2.set_ylabel(r'cross section [mb]')
y_label2 = (r"$\frac{3}{2 \pi} \lambda_0^2 \left( \frac{\Gamma}{\Gamma + \Gamma'} \right)$", r"$\frac{3}{2 \pi} \lambda_0^2 \left( \frac{\Gamma}{\Gamma + \Gamma'} \right)^2$", r"$\sigma_{Th}$")
y_tick = [barn*3/(2*np.pi)*lamda0**2*L/(L+L1), barn*3/(2*np.pi)*lamda0**2*(L/(L+L1))**2, barn*sigma_th()]
plt.yticks(y_tick, y_label2)

ax3 = ax1.secondary_xaxis('bottom')
# show result
pyl.show()


