import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22}) # larger font size
from math import pi

# value of first m terms at x,t
# with respect to sum this is k=0,1,...,m-1
def SeriesTerms(x,t,m):
    Aux = 0
    for k in range(m):
        Aux += 8*(-1)**k/(pi*(2*k+1)**2) * np.exp(-(2*k+1)**2*t/4) * np.sin((2*k+1)*x/2)
    return Aux
    
Nx = 1000
xvals = np.linspace(0,pi,Nx)
mvals = [1,2,5,10,100]

fig,ax = plt.subplots(1,3,figsize=(18,6))
# time 0
max_error = []
for m in mvals:
    yvals = [SeriesTerms(x,0,m) for x in xvals]
    ax[0].plot(xvals,yvals)
    max_error.append(max([abs(yvals[i]-xvals[i]) for i in range(Nx)]))
ax[0].legend(['m='+str(m) for m in mvals])
ax[0].set_title(r't=0')
ax[0].set_xlabel('x')
ax[0].set_ylabel('u(x,0)')

# show maximum error
ax[2].plot(mvals,max_error)
ax[2].set_title('max initial error')
ax[2].set_xlabel('m')
ax[2].set_ylabel('max error')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
print(max_error[-1])

# time 2
for m in mvals:
    yvals = [SeriesTerms(x,2,m) for x in xvals]
    ax[1].plot(xvals,yvals)
ax[1].set_title(r't=2')
ax[1].set_xlabel('x')
ax[1].set_ylabel('u(x,2)')


plt.tight_layout()
plt.savefig('Heat_Series.pdf')