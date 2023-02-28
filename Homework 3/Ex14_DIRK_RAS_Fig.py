# Compute the region of absolute stability for DIRK2
# adapted from ODEsolver notes
import numpy as np
import matplotlib.pyplot as plt


# stability function
gamma  = 1-1/(2**.5)
R = lambda z: (1+(1-2*gamma)*z)/((1-gamma*z)*(1-gamma*z))

# Construct a mesh
nx = 500
ny = 500
x = np.linspace(-5,5,nx)
y = np.linspace(-5,5,ny)
xg,yg = np.meshgrid(x,y)

# Evaluate
z = xg + 1j*yg
f = R(z)
absf = (((f.real)**2 + (f.imag)**2) < 1)*2 # set unstable region to 2, stable to 0

# plot
plt.rcParams.update({'font.size': 22})
fig, ax = plt.subplots(figsize=(8,8))
plt.contourf(xg,yg,absf,np.arange(2))
plt.title("RAS")
plt.xlabel("Re(z)")
plt.ylabel("Im(z)")
ax.set_aspect(1)
plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.savefig('DIRK_P2_RAS.pdf')