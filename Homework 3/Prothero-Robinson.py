#Packages
from math import sin,cos,pi,exp
import numpy as np
import matplotlib.pyplot as plt

### ODE info ###
# params
phi  = lambda t: sin(t+pi/4)
Dphi = lambda t: cos(t+pi/4)
L    = 10000
Tmax = 10
# Exact soln
y = lambda t: exp(-L*t)*(y0 - phi(0)) + phi(t)

### DIRK solvers ###
gamma = 1 - 1/(2**.5)
def DIRK2_Step(tn,un,h):
    k1 = (-L*(un-phi(tn+gamma*h))+Dphi(tn+gamma*h))/(1+gamma*h*L)
    k2 = (-L*(un+(1-gamma)*h*k1-phi(tn+h))+Dphi(tn+h))/(1+gamma*h*L)
    un_new = un + h*(1-gamma)*k1 + h*gamma*k2
    return un_new
    
def DIRK3_Step(tn,un,h):
    k1 = (-L*(un-phi(tn+gamma*h))+Dphi(tn+gamma*h))/(1+gamma*h*L)
    k2 = (-L*(un+(1-2*gamma)*h*k1-phi(tn+(1-gamma)*h))+Dphi(tn+(1-gamma)*h))/(1+gamma*h*L)
    un_new = un+h*k1/2+h*k2/2
    return un_new


### Part a compute each solution ###
y0 = sin(pi/4)
h_vec = 10**(-np.arange(1,6,5/24)) # time-step array
error2 = [] # DIRK2 error
error3 = [] # DIRK3 error
for h in h_vec:
    Times = np.arange(0,Tmax,h)
    Y_vec = [y(t) for t in Times] # exact soln.
    Y_approx2 = [y0]
    Y_approx3 = [y0]
    for i in range(len(Times)-1):
        Y_approx2.append(DIRK2_Step(Times[i],Y_approx2[-1],h))
        Y_approx3.append(DIRK3_Step(Times[i],Y_approx3[-1],h))
    error2.append(max([abs(Y_approx2[i]-Y_vec[i]) for i in range(len(Times))]))
    error3.append(max([abs(Y_approx3[i]-Y_vec[i]) for i in range(len(Times))]))
    
# plot
plt.rcParams.update({'font.size': 22})
fig, ax = plt.subplots(figsize=(8,8))
plt.plot(h_vec, [.00005*h for h in h_vec],linestyle = '-',color='grey', linewidth=3)
plt.plot(h_vec, [.046*h**2 for h in h_vec],linestyle = '--',color='grey', linewidth=3)
plt.plot(h_vec,error2, linewidth=3)
plt.plot(h_vec,error3, linewidth=3)
plt.xlabel("step-size (h)")
plt.ylabel("error (e)")
plt.xscale('log')
plt.yscale('log')
plt.legend(['slope 1','slope 2','DIRK2','DIRK3'])
plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('Proth_Rob_a.pdf')


### part b ###
y0 = sin(pi/4)+10
fig, ax = plt.subplots(nrows=2,ncols=2,figsize=(16,16))
for h in [10**-1,10**-2,10**-3]:
    for (row,Tmax) in enumerate([10,1]):
        # Tmax = 10
        Times = np.arange(0,Tmax,h)
        Y_vec = np.array([y(t) for t in Times]) # exact soln.
        Y_approx2 = [y0]
        Y_approx3 = [y0]
        for i in range(len(Times)-1):
            Y_approx2.append(DIRK2_Step(Times[i],Y_approx2[-1],h))
            Y_approx3.append(DIRK3_Step(Times[i],Y_approx3[-1],h))
        Y_approx2 = np.array(Y_approx2)
        error2 = np.abs(Y_vec-Y_approx2)
        ax[row,0].plot(Times[1:],error2[1:])
        Y_approx3 = np.array(Y_approx3)
        error3 = np.abs(Y_vec-Y_approx3)
        ax[row,1].plot(Times[1:],error3[1:])
    
# plot
ax[0,0].set_title('DIRK2, Tmax=10')
ax[0,1].set_title('DIRK3, Tmax=10')
ax[1,0].set_title('DIRK2, Tmax=1')
ax[1,1].set_title('DIRK3, Tmax=1')
ax[1,1].legend([r'$h=10^{-1}$',r'$h=10^{-2}$',r'$h=10^{-3}$'])
for (i,j) in [(0,0),(0,1),(1,0),(1,1)]:
    ax[i,j].set_yscale('log')
    ax[i,0].set_ylabel('Error |e(t)|')
    ax[1,j].set_xlabel('Time (t)')
    ax[i,j].grid(color='k', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('Proth_Rob_b.pdf')