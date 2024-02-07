import numpy as np
from numba import njit

@njit
def standard_map(q, p, K=0.6):
    pnew = (p + K * np.sin(q))
    qnew = (q + pnew)
    return qnew, pnew
