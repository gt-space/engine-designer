# --------------------------------------------
#        TRANSONIC ANALYSIS (TA) MODULE
# --------------------------------------------


import numpy as np
import scipy as sp

N = 50 #number of zones, more zones --> more accurate simulation

# I don't know what the Mdots are rn, Anushka's getting it out, but it should be a vector of length n-1
Mdot = []

mDot = np.zeros(N)
rhoStar = np.zeros(N)
aStar = np.zeros(N)

Yn = zoneY(N, mDot, rhoStar, aStar)

def zoneY(N, mDot, rhoStar, aStar):
    # also don't know K, A but we can calculate them with some variables Anushka is getting
    K = []
    A = []
    
    for i in range(N-1):
        K(i) = mDot(i+1)/mDot(i)
        A(i) = rhoStar(i) * aStar(i) / (rhoStar(i+1) * aStar(i+1))

    tridiagArr = np.zeros(N-1,N-1)
    for i in range(N-1):
        tridiagArr[i,i] = 1 + A[i] * K[i]
        if i > 0:
            tridiagArr[i,i-1] = -A[i] * K[i]
        if i < N-2:
            tridiagArr[i,i+1] = -1

    bArr = np.zeros(N-1)
    bArr(-1) = 1

    Yn = np.linalg.solve(tridiagArr, bArr) ** 0.5
    return Yn

def initLine(N, rw, zw, zaxis):
    # zw is 
    # zaxis is found such that the mach number at r = 0, z = zaxis is the same as at the point r = rw, z = zw
    r = (rw * np.sin(1/N * np.pi/2))^1.2
    z = zw + ((rw - r)/ rw)^2 * (zaxis - zw)
    return r, z

