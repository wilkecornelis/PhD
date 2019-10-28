import numpy as np


shape = (10,1)
x = np.ones(shape)
y = np.ones(shape)

g = 5*np.ones(shape)
g[4] = 45

a = g*np.exp(3*x + 4*y)

print(a)