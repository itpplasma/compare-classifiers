import numpy as np
import matplotlib.pyplot as plt

TOL_TIPS = 1e-12

# Data
z_flux = np.loadtxt('orbit.dat')
z_can = np.loadtxt('orbit_can.dat')

var_cut = np.loadtxt('cut.dat')
var_cut_can = np.loadtxt('cut_can.dat')

index_tips = np.abs(var_cut[:, 4]) < TOL_TIPS
var_tips = var_cut[index_tips, :]
var_tips_can = var_cut_can[index_tips, :]

index_torcut = np.logical_not(index_tips)
var_torcut = var_cut[index_torcut, :]
var_torcut_can = var_cut_can[index_torcut, :]

plt.figure()
plt.plot(var_tips[:, 0], var_tips[:, 1], 'r.')
plt.xlabel('s')
plt.ylabel('theta')

plt.figure()
plt.plot(var_torcut[:, 0], var_torcut[:, 1], 'b.')
plt.xlabel('s')
plt.ylabel('theta')

# Canonical
plt.figure()
plt.title('Canonical variables')
plt.plot(var_tips_can[:, 1], var_tips_can[:, 0], 'r.')
plt.xlabel('p_theta')
plt.ylabel('theta')

plt.figure()
plt.title('Canonical variables')
plt.plot(var_torcut_can[:, 1], var_torcut_can[:, 0], 'b.')
plt.xlabel('p_theta')
plt.ylabel('theta')
