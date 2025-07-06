import numpy as np
import matplotlib.pyplot as plt


fig, ax = plt.subplots(1, 1, figsize = (5, 5))
N = 100
y = np.full(N, 0.)
x = np.linspace(0, 1.0, N)

ax.plot(x, y)