import numpy as np
import matplotlib.pyplot as plt

# Load the standard_map function (assuming it's in a separate file named standard_map.py)
from standard_map import standard_map

# Number of iterations and initial conditions
tmax = 1000
norb = 100

# Arrays to hold the results
q = np.empty((norb, tmax + 1))
p = np.empty((norb, tmax + 1))

# Initialize q to cut all orbits and p on a regular grid
q[:, 0] = np.pi
p[:, 0] = np.linspace(0, 2 * np.pi, num=norb)

for k in range(tmax):
    for i in range(norb):
        q[i, k + 1], p[i, k + 1] = standard_map(q[i, k], p[i, k])

# Plot the results
plt.scatter(q, p, s=0.1)
plt.xlabel('q')
plt.ylabel('p')
plt.title('Standard Map')
plt.show()
