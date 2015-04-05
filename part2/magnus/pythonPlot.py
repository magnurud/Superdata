import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import linalg, optimize
from sympy.mpmath import plot
from mpl_toolkits.mplot3d import Axes3D

with open('meh.asc') as f:
	a = f.readline().split()
with open('X.asc') as f:
	b = f.readline().split()
with open('Y.asc') as f:
	c = f.readline().split()
        
print float(c[3])
U = []
X = []
Y = []
for i in range(len(a)):
    U.append(float(a[i]))
    X.append(float(b[i]))
    Y.append(float(c[i]))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, U)

plt.show()

