import numpy as np
from numba import njit

@njit
def standard_map(q, p, K=0.6):
    pnew = np.mod(p + K * np.sin(q), 2 * np.pi)
    qnew = np.mod(q + pnew, 2 * np.pi)
    return qnew, pnew
