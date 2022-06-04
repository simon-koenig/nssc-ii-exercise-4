import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

# http://en.wikipedia.org/wiki/File:Bilininterp.png
xi = np.array([0.0, 1.0])
yi = np.array([0.0, 1.0])
zi = np.array([[0.0, 1.0], [1.0, 0.5]])

# Another example
xi = np.array([0.0, 0.25, 1.0])
yi = np.array([0.0, 0.75, 1.0])
zi = np.array([[0.0, 0.5, 1.0], [0.5, 0.7, 0.5], [1.0, 1.0, 1.0]])

# I want 20 "levels" to be shown
contour_breaks = 20
ticks = np.linspace(zi.min(), zi.max(), contour_breaks, endpoint=True)

# Attempt 4 (interp2d does to correct bilinear interpolation)
fig = plt.figure()
axes = fig.add_subplot(111, aspect='equal')
f = interp2d(xi, yi, zi, kind='linear')
xi2 = np.linspace(0., 1., 100)
yi2 = np.linspace(0., 1., 100)
zi2 = f(xi2, yi2)
axes.contour(xi2, yi2, zi2, ticks[1:-1], colors='k')
fill = axes.contourf(xi2, yi2, zi2, ticks, cmap=cm.jet)
fig.colorbar(fill, ticks=ticks)

# Show the plots
plt.show()