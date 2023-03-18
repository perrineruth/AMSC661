# packages
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt

# build matrices
main  = [20]*3 + [11] + [2]*5
lower = [-10]*3 + [-1]*5
upper = [-10]*3 + [-1]*5
# sparse array for good habits
A = sparse.diags([lower,main,upper],[-1,0,1]).tocsc()
b = [0]*8+[10]
u = sparse.linalg.spsolve(A,b)
# add begin and end
u = np.concatenate(([0], u, [10]))
# print(u)

Actual = [5/32*n for n in range(4)] + [10 + 25/16*(n-10) for n in range(4,10+1)]
# print(Actual)

# plot
plt.rcParams.update({'font.size': 22})
fig, ax = plt.subplots(figsize=(8,8))
plt.plot(range(10+1),u,'x')
plt.plot(range(10+1),Actual)
plt.xlabel("x")
plt.ylabel("u")
plt.legend(('numerical','analytic'))
ax.set_aspect(1)
plt.grid(color='k', linestyle='--', linewidth=0.5)
plt.savefig('HeatEqP1.pdf')