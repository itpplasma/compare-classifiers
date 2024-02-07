import numpy as np
from scipy.optimize import minimize
import random
import scipy
from standard_map import standard_map
import matplotlib.pyplot as plt
from symp_gpr import sympGPR

# init parameters
K = 0.6      # stochasticity parameter
N = 20        # training data
nm = 100      # map applications
Ntest = 50    # test data
sig2_n = 1e-12 #noise**2 in observations

# sample initial data points from halton sequence in (0, 2\pi) x (0, 2 \pi)
seq = scipy.stats.qmc.Halton(2)
X0 = seq.random(N)*np.array([2*np.pi, 2*np.pi])

# Arrays to hold the results for training data
q = np.empty((N, 2))
p = np.empty((N, 2))

# Initialize q and p
q[:, 0] = X0[:, 0]
p[:, 0] = X0[:, 1]

# calculate one iteration of standard map to get training data
for k in range(1):
    for i in range(N):
        q[i, k + 1], p[i, k + 1] = standard_map(q[i, k], p[i, k], K)

#set training data
class Set_data:
    def __init__(self, q, p):
        self.q = q[:, 0]
        self.p = p[:, 0]
        self.Q = q[:, 1]
        self.P = p[:, 1]

# Set test data on regular grid
qtest = np.empty((Ntest, nm + 1))
ptest = np.empty((Ntest, nm + 1))

# Initialize q to cut all orbits and p on a regular grid
qtest[:, 0] = np.pi
ptest[:, 0] = np.linspace(0, 2 * np.pi, num=Ntest)

for k in range(nm):
    for i in range(Ntest):
        qtest[i, k + 1], ptest[i, k + 1] = standard_map(qtest[i, k], ptest[i, k], K)
X0test = np.stack((np.mod(qtest, 2*np.pi), np.mod(ptest-np.pi, 2*np.pi)))

trainingdata = Set_data(q, p)
testdata = Set_data(qtest, ptest)

# plot training data and test data
plt.figure()
plt.plot(q[:, 0], p[:, 0], 'rx')
plt.plot(q[:, 1], p[:, 1], 'kx')
plt.title('Training data')

# call SympGPR 
method = 'explicit'

outq, outp = sympGPR(method, trainingdata, testdata, sig2_n, nm)

pmap = np.mod(outp-np.pi, 2*np.pi)
qmap = np.mod(outq, 2*np.pi)
fmap = np.stack((qmap, pmap))

plt.figure(figsize = [10,3])
plt.subplot(1,3,1)
for i in range(0, Ntest):
    plt.plot(qmap[:,i], pmap[:,i], 'k^', label = 'GP', markersize = 0.5)
plt.xlabel(r"$\theta$", fontsize = 20)
plt.ylabel(r"I", fontsize = 20)
plt.tight_layout()

plt.subplot(1,3,2)
plt.plot(X0test[0], X0test[1],  color = 'dodgerblue', marker = 'o', linestyle = 'None',  markersize = 0.5)
plt.xlabel(r"$\theta$", fontsize = 20)
plt.ylabel(r"I", fontsize = 20)
plt.tight_layout()

plt.subplot(1,3,3)
plt.plot(X0test[0], X0test[1],  color = 'dodgerblue', marker = 'o', linestyle = 'None',  markersize = 0.5)
for i in range(0, Ntest):
    plt.plot(qmap[:,i], pmap[:,i], 'k^', label = 'GP', markersize = 0.5)
plt.xlabel(r"$\theta$", fontsize = 20)
plt.ylabel(r"I", fontsize = 20)
plt.tight_layout()