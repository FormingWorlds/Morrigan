import numpy as np 
from pylaplace import LaplaceCoefficient 
from constants import *
import pdb

def secular_solution(ap, Mp, ecc, Rp, Ms, N):
    varpi = np.random.uniform(0.0, 2.0 * np.pi, N)

    h0 = ecc * np.sin(varpi)
    k0 = ecc * np.cos(varpi)

    mean_motion = np.sqrt(G*Ms/ap**3)
    A = np.zeros((N,N)) #empty interaction matrix
    apo = (1.0 + ecc) * ap
    peri = (1.0 - ecc) * ap

    laplace = LaplaceCoefficient(method = 'Brute') #to calculate laplace coefficients
    for i in range(N): 
        for j in range(N): 
            if i == j:
                continue #skip self-interactions
            if min(apo[i], apo[j]) > max(peri[i], peri[j]):
                continue #skip overlapping orbits
            if ap[i] < ap[j]:
                alpha = ap[i]/ap[j]
                alpha_bar = alpha 
            else:
                alpha = ap[j]/ap[i]
                alpha_bar = 1

            #calculate laplacian coefficients, b1,3/2 and b2,3/2
            #(a,s,m,p,q)
            coeff_m1 = laplace(alpha, 3/2, 1, 1, 1) #m = 1 for A_ii
            coeff_m2 = laplace(alpha, 3/2, 2, 1, 1) #m = 2 for A_ij

            factor = mean_motion[i] * 0.25 * Mp[j]/(Ms + Mp[i]) * alpha * alpha_bar
            A[i,i] += factor * coeff_m1
            A[i,j] = -factor * coeff_m2

    #solve the eigenvalue problem
    pdb.set_trace()
    g, S = np.linalg.eig(A)
    g = np.real(g)
    S = np.real(S)

    #solve for integration constants using S
    Csinb = np.linalg.solve(S, h0)
    Ccosb = np.linalg.solve(S, k0)

    #calculate amplitudes (C) and phase angles (beta)
    C = np.sqrt(Csinb**2 + Ccosb**2)
    beta = np.arctan2(Csinb, Ccosb)

    #scale eigenvectors by the amplitudes (columns of ecc_vec)
    ecc_vec = S * C

    return ecc_vec, g, beta