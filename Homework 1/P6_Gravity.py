import numpy as np
from math import pi, cos, sin
import matplotlib.pyplot as plt

T = 4*pi

# compute derivative according to gravity equation
def D_func(prev):
    return np.array([prev[2],\
                     prev[3],\
                     -prev[0]/(prev[0]**2+prev[1]**2),\
                     -prev[1]/(prev[0]**2+prev[1]**2)])


# plot exact solution
Times = [2*pi/100 * i for i in range(101)]
plt.plot(np.cos(Times),np.sin(Times))

# number of iterations instead of
for N in [40,80,160]:
    # time step
    h = T/N
    # all time values
    Times  = [i * h for i in range(N+1)]
    
    # state variable vector
    States = np.zeros((4,N+1))
    
    # initial conditions
    States[:,0] = np.array([1,0,0,1])
    # second state exact
    States[:,1] = np.array([cos(h),sin(h),-sin(h),cos(h)])
    
    # iterate through time
    for t in range(2,5):
        States[:,t] = -4*States[:,t-1] + 5*States[:,t-2] + h*(4*D_func(States[:,t-1]) + 2*D_func(States[:,t-2]))
        
    print(np.linalg.norm(States[:,0] - States[:,-1]))
    # plot simulation
    plt.plot(States[0,:],States[1,:])
    plt.xlim([ -2, 2])
    plt.ylim([ -2, 2])
    
    
plt.legend(['Exact','N=20','N=40','N=80'])
plt.show()