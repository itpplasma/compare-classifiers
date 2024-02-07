import numpy as np
from scipy.optimize import minimize
import scipy
from scipy.optimize import newton, bisect
from scipy.linalg import solve_triangular

# from kernels import * # implicit SympGPR
from kernels_expl_per_q_sq_p import * #explicit SympGPR
def f_kern(x, y, x0, y0, l):
    return kern_num(x,y,x0,y0,l[0], l[1])

def d2kdxdx0(x, y, x0, y0, l):
    return d2kdxdx0_num(x,y,x0,y0,l[0], l[1])

def d2kdydy0(x, y, x0, y0, l):
    return d2kdydy0_num(x,y,x0,y0,l[0], l[1])

def d2kdxdy0(x, y, x0, y0, l):
    return d2kdxdy0_num(x,y,x0,y0,l[0], l[1])

def d2kdydx0(x, y, x0, y0, l):
    return d2kdxdy0(x, y, x0, y0, l)


def nll_transform(log10hyp, sig, sig2n, x, y, N):
    hyp = 10**log10hyp
    return nll_chol(np.hstack((hyp, sig, [sig2n])), x, y, N)

def guessP(x, y, hypp, xtrainp, ztrainp, Kyinvp, N):
    Ntest = 1
    Kstar = np.empty((Ntest, int(len(xtrainp)/2)))
    buildKreg(np.hstack((x,y)), xtrainp, hypp, Kstar)
    Ef = Kstar.dot(Kyinvp.dot(ztrainp))
    return Ef

def calcQ(x,y, xtrain, l, Kyinv, ztrain):
    # get \Delta q from GP on mixed grid.
    Kstar = np.empty((len(xtrain), 2))
    build_K(xtrain, np.hstack(([x], [y])), l, Kstar)
    qGP = Kstar.T.dot(Kyinv.dot(ztrain))
    dq = qGP[1]
    return dq

def Pnewton(P, x, y, l, xtrain, Kyinv, ztrain):
    Kstar = np.empty((len(xtrain), 2))
    build_K(xtrain, np.hstack((x, P)), l, Kstar)
    pGP = Kstar.T.dot(Kyinv.dot(ztrain))
    f = pGP[0] - y + P
    # print(pGP[0])
    return f

def calcP(x,y, l, hypp, xtrainp, ztrainp, Kyinvp, xtrain, ztrain, Kyinv, Ntest):
    # as P is given in an implicit relation, use newton to solve for P (Eq. (42))
    # use the GP on regular grid (q,p) for a first guess for P
    pgss = guessP([x], [y], hypp, xtrainp, ztrainp, Kyinvp, Ntest)
    res, r = newton(Pnewton, pgss, full_output=True, maxiter=205000, disp=True,
        args = (np.array([x]), np.array ([y]), l, xtrain, Kyinv, ztrain))
    return res

def calcP_expl(x,y, l, xtrain, ztrain, Kyinv):
    Kstar = np.empty((len(xtrain), 2))
    build_K(xtrain, np.hstack((x,y)), l, Kstar)
    pGP = Kstar.T.dot(Kyinv.dot(ztrain))
    res = -pGP[0] + y
    return res

def buildKreg(xin, x0in, hyp, K):
    # set up "usual" covariance matrix for GP on regular grid (q,p)
    l = hyp[:-1]
    sig = hyp[-1]
    N = K.shape[0]
    N0 = K.shape[1]
    x0 = x0in[0:N0]
    x = xin[0:N]
    y0 = x0in[N0:2*N0]
    y = xin[N:2*N]
    for k in range(N):
        for lk in range(N0):
            K[k,lk] = f_kern(
                     x0[lk], y0[lk], x[k], y[k], l)
    K[:,:] = sig*K[:,:]

def build_K(xin, x0in, hyp, K):
    # set up covariance matrix with derivative observations, Eq. (38)
    l = hyp[:-1]
    sig = hyp[-1]
    N = K.shape[0]//2
    N0 = K.shape[1]//2
    x0 = x0in[0:N0]
    x = xin[0:N]
    y0 = x0in[N0:2*N0]
    y = xin[N:2*N]
    for k in range(N):
        for lk in range(N0):
            K[k,lk] = d2kdxdx0(
                x0[lk], y0[lk], x[k], y[k], l)
            K[N+k,lk] = d2kdxdy0(
                 x0[lk], y0[lk], x[k], y[k], l)
            K[k,N0+lk] = d2kdydx0(
                 x0[lk], y0[lk], x[k], y[k], l)
            K[N+k,N0+lk] = d2kdydy0(
                x0[lk], y0[lk], x[k], y[k], l)
    K[:,:] = sig*K[:,:]

def nll_chol(hyp, x, y, N):
    K = np.empty((N, N))
    build_K(x, x, hyp[:-1], K)
    Ky = K + np.abs(hyp[-1])*np.diag(np.ones(N))
    L = scipy.linalg.cholesky(Ky, lower = True)
    alpha = solve_cholesky(L, y)
    ret = 0.5*y.T.dot(alpha) + np.sum(np.log(L.diagonal()))
    return ret

def build_K_expl(xin, x0in, hyp, K):
    # set up covariance matrix with derivative observations, Eq. (38)
    l = hyp[:-1]
    sig = hyp[-1]
    N = K.shape[0]//2
    N0 = K.shape[1]//2
    x0 = x0in[0:N0]
    x = xin[0:N]
    y0 = x0in[N0:2*N0]
    y = xin[N:2*N]
    for k in range(N):
        for lk in range(N0):
            K[k,lk] = d2kdxdx0(
                x0[lk], y0[lk], x[k], y[k], l)
            K[N+k,lk] = d2kdxdy0(
                 x0[lk], y0[lk], x[k], y[k], l)
            K[k,N0+lk] = d2kdydx0(
                 x0[lk], y0[lk], x[k], y[k], l)
            K[N+k,N0+lk] = d2kdydy0(
                x0[lk], y0[lk], x[k], y[k], l)
    K[:,:] = sig*K[:,:]

def nll_expl(hyp, x, y, N, ind):
    K = np.empty((N, N))
    if ind == 0:
        build_K_expl(x, x, np.hstack((hyp[0], 0, hyp[1])), K)
        Ky = K + np.abs(hyp[-1])*np.diag(np.ones(N))
        Ky = Ky[0:len(y), 0:len(y)]
    else:
        build_K_expl(x, x, np.hstack((0, hyp[0], hyp[1])), K)
        Ky = K + np.abs(hyp[-1])*np.diag(np.ones(N))
        Ky = Ky[len(y):2*len(y), len(y):2*len(y)]

    L = scipy.linalg.cholesky(Ky, lower = True)
    alpha = solve_cholesky(L, y)
    nlp_val = 0.5*y.T.dot(alpha) + np.sum(np.log(L.diagonal()))

    return nlp_val

def solve_cholesky(L, b):
    return solve_triangular(
        L.T, solve_triangular(L, b, lower=True, check_finite=False),
        lower=False, check_finite=False)

def applymap(nm, Ntest, l, hypp, Q0map, P0map, xtrainp, ztrainp, Kyinvp, xtrain, ztrain, Kyinv):
    # Application of symplectic map
    #init
    pmap = np.zeros([nm, Ntest])
    pdiff = np.zeros([nm, Ntest])
    qmap = np.zeros([nm, Ntest])
    #set initial conditions
    pmap[0,:] = P0map
    qmap[0,:] = Q0map
    pdiff[0,:] = P0map
    # loop through all test points and all time steps
    for i in range(0,nm-1):
        for k in range(0, Ntest):
            # set new P including Newton for implicit Eq. (42)
            pmap[i+1, k] = calcP(qmap[i,k], pmap[i, k], l, hypp, xtrainp, ztrainp, Kyinvp, xtrain, ztrain, Kyinv, Ntest)

            pdiff[i+1, k] = pdiff[i, k] + (pmap[i+1, k] - pmap[i, k])
            pmap[i+1, k] = np.mod(pmap[i+1, k], 2*np.pi) # only for standard map

        for k in range(0, Ntest):
            if np.isnan(pmap[i+1, k]):
                qmap[i+1,k] = np.nan
            else:
                # then: set new Q via calculating \Delta q and adding q (Eq. (43))
                dqmap = calcQ(qmap[i,k], pmap[i+1,k], xtrain, l, Kyinv, ztrain)
                qmap[i+1, k] = np.mod(dqmap + qmap[i, k], 2.0*np.pi)
    return qmap, pmap, pdiff

def nll_transform_expl(log10hyp, sig, sig2n, x, y, N, ind):
    hyp = 10**log10hyp
    # hyp = log10hyp
    # print(hyp)
    out = nll_expl(np.hstack((hyp, sig, [sig2n])), x, y, N, ind)
    return out

def applymap_expl(nm, Ntest, l, Q0map, P0map, xtrain, ztrain, Kyinv):
    # Application of symplectic map
    #init
    pmap = np.zeros([nm, Ntest])
    pdiff = np.zeros([nm, Ntest])
    qmap = np.zeros([nm, Ntest])
    #set initial conditions
    pmap[0,:] = P0map
    qmap[0,:] = Q0map
    pdiff[0,:] = P0map

    # loop through all test points and all time steps
    for i in range(0,nm-1):
        for k in range(0, Ntest):
            # set new P including Newton for implicit Eq (42)
            # print('nm = ', i, 'N_i = ',k)
            pmap[i+1, k] = calcP_expl(qmap[i,k], pmap[i, k], l, xtrain, ztrain, Kyinv)
            pdiff[i+1, k] = pdiff[i, k] + (pmap[i+1, k] - pmap[i, k])
            pmap[i+1, k] = np.mod(pmap[i+1, k], 2*np.pi)
        for k in range(0, Ntest):
            if np.isnan(pmap[i+1, k]):
                qmap[i+1,k] = np.nan
            else:
                # then: set new Q via calculating \Delta q and adding q (Eq. (43))
                dqmap = calcQ(qmap[i,k], pmap[i+1,k], xtrain, l, Kyinv, ztrain)
                qmap[i+1, k] = dqmap + qmap[i, k]
    return qmap, pmap, pdiff
    
def sympGPR(method, trainingdata, testdata, sig2_n, nm):
    # set up GP
    # as indicated in Algorithm 1: Semi-implicit symplectic GP map
    # hyperparameter optimization of length scales (lq, lp)
    q = trainingdata.q
    p = trainingdata.p
    Q = trainingdata.Q
    P = trainingdata.P
    q0test = testdata.q
    p0test = testdata.p

    N = len(q)
    print('N =', N)
    Ntest = len(q0test)

    zqtrain = Q - q
    zptrain = p - P

    xtrain = q.flatten()
    ytrain = P.flatten()
    xtrain = np.hstack((q, P)).T
    ztrain = np.concatenate((zptrain.flatten(), zqtrain.flatten()))

    if method == 'implicit':
        log10l0 = np.array((-1, -1), dtype = float)
        
        #  Step 1: Usual GP regression of P over (q,p)
        #fit GP + hyperparameter optimization to have a first guess for newton for P
        xtrainp = np.hstack((q, p))
        ztrainp = P - p
        
        sigp = 2*np.amax(np.abs(ztrainp))**2
        def nll_transform2(log10hyp, sig, sig2n, x, y, N):
            hyp = 10**log10hyp
            return nll_chol(np.hstack((hyp, sig, [sig2n])), x, y, N)
    
        res = minimize(nll_transform2, np.array((log10l0)), args = (sigp, sig2_n, xtrainp, ztrainp.T.flatten(), N), method='L-BFGS-B', bounds = ((-10, 1), (-10, 1)))
        
        lp = 10**res.x
        hypp = np.hstack((lp, sigp))
        print('Optimized lengthscales for regular GP: lq =', "{:.2f}".format(lp[0]), 'lp = ', "{:.2f}".format(lp[1]))
        
        # build K and its inverse
        Kp = np.zeros((N, N))
        buildKreg(xtrainp, xtrainp, hypp, Kp)
        Kyinvp = scipy.linalg.inv(Kp + sig2_n*np.eye(Kp.shape[0]))
        #%%
        # Step 2: symplectic GP regression of -Delta p and Delta q over mixed variables (q,P) according to Eq. 41
        # hyperparameter optimization for lengthscales (lq, lp) and GP fitting
        sig = 2*np.amax(np.abs(ztrain))**2
        log10l0 = np.array((0,0), dtype = float)
        
        res = minimize(nll_transform, np.array((0,-1)), args = (sig, sig2_n, xtrain, ztrain.T.flatten(), 2*N), method='L-BFGS-B', tol= 1e-8, bounds = ((-2, 2), (-2, 2)))#, 
        
        sol1 = 10**res.x
        
        l = [np.abs(sol1[0]), np.abs(sol1[1])]
        print('Optimized lengthscales for mixed GP: lq =', "{:.2f}".format(l[0]), 'lp = ', "{:.2f}".format(l[1]), 'sig = ', "{:.2f}".format(sig))
        
        #build K(x,x') and regularized inverse with sig2_n
        # K(x,x') corresponds to L(q,P,q',P') given in Eq. (38) 
        hyp = np.hstack((l, sig))
        K = np.empty((2*N, 2*N))
        build_K(xtrain, xtrain, hyp, K)
        Kyinv = scipy.linalg.inv(K + sig2_n*np.eye(K.shape[0]))
        
        # caluate training error
        # Eftrain = K.dot(Kyinv.dot(ztrain))
        # outtrain = mean_squared_error(ztrain, Eftrain)
        # print('training error', "{:.1e}".format(outtrain))
        
        #% Application of symplectic map
        outq, outp, pdiff = applymap(
            nm, Ntest, hyp, hypp, q0test, p0test, xtrainp, ztrainp, Kyinvp, xtrain, ztrain, Kyinv)
        
    #%
    elif method == 'explicit':
        # use kernels_expl_per_q_sq_p.pyd for sum kernel in func.py

        # Step 2: symplectic GP regression of -Delta p and Delta q over mixed variables (q,P) according to Eq. 41
        # hyperparameter optimization for lengthscales (lq, lp) and GP fitting
        # this is the explicit method
        # lq and lp are trained separately

        sig = 2*np.amax(np.abs(ztrain))**2
        log10l0 = np.array((1), dtype = float)
       
        #log 10 -> BFGS
        res_lq = minimize(nll_transform_expl, np.array((1)), args = (sig, 1e-8, xtrain, zptrain.T.flatten(), 2*N, 0), method='L-BFGS-B')#, bounds = (-2, 2))#, 
        res_lp = minimize(nll_transform_expl, np.array((1)), args = (sig, 1e-8, xtrain, zqtrain.T.flatten(), 2*N, 1), method='L-BFGS-B')#, bounds = ((-2, 2), (-2, 2)),tol= 1e-8)#, 

        sol1 = 10**res_lq.x
        sol2 = 10**res_lp.x

        l = np.hstack((np.abs(sol1), np.abs(sol2)))
        print('Optimized lengthscales for mixed GP: lq =', "{:.2f}".format(l[0]), 'lp = ', "{:.2f}".format(l[1]))
        hyp = np.hstack((l, sig))
        K = np.empty((2*N, 2*N))
        build_K(xtrain, xtrain, hyp, K)
        Kyinv = scipy.linalg.inv(K + sig2_n*np.eye(K.shape[0]))
        
        # caluate training error
        # Eftrain = K.dot(Kyinv.dot(ztrain))
        # outtrain = mean_squared_error(ztrain, Eftrain)
        # print('training error', "{:.1e}".format(outtrain))
        
        outq, outp, pdiff = applymap_expl(
        nm, Ntest, hyp, q0test, p0test, xtrain, ztrain, Kyinv)

    return outq, outp