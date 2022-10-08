# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 13:53:53 2022

@author: coolc
"""

import numpy as np


n = 3000
l = 10

A = np.eye(n) + np.random.rand(n,n) * 1e-5

V = np.zeros((n, l))

for i in range(l):
    V[:, i] = np.random.rand(n)


Tol = 1e-3
iter = 0
while(True)    :
    S = V.T @ A @ V
    
    lamb, z = np.linalg.eig(S)[0], np.linalg.eig(S)[1]
    
    
    x = np.zeros((n, l))
    
    for k in range(l):
        sum = None
        for i in range(l):
            if sum is None:
                sum = z[i, k] * V[:, i]
            else:
                sum += z[i, k] * V[:, i]
        x[:, k] = sum
        
        
    # calculating residual vector
    r = np.zeros((n, l))
    I = np.eye(n)
    
    for k in range(l):
        r[:, k] = (A - lamb[k] * I) @ x[:, k]
    
    
    # calculate correction vector
    
    d = np.zeros((n, l))
    A_ii = np.zeros_like(A)
    
    for i in range(A.shape[0]):
        A_ii[i, i] = A[i, i] 
    
    for k in range(l):
        d[:, k] = np.linalg.inv(lamb[k] - A_ii) @ r[:, k]
        
    # Expand the subspace C(V) with the orthonormalized correction vectors. 
    # First, project on the orthogonal complement of C(V)
    
    q = (np.eye(n) - V @ V.T) @ d    
    
    q_norm = np.linalg.norm(q)
    
    print(f"iter {iter} : {q_norm}")
    
    v_p = q / q_norm
    
    
    V = v_p
    
    if q_norm < Tol:
        break
    if iter > 100:
        break
    iter += 1

# # #  Davidson's method : 다른 사람 예제

# from __future__ import division
# from __future__ import print_function
# import math
# import numpy as np
# import time

# ''' Block Davidson, Joshua Goings (2013)

#     Block Davidson method for finding the first few
#  	lowest eigenvalues of a large, diagonally dominant,
#     sparse Hermitian matrix (e.g. Hamiltonian)
# '''

# n = 3000					# Dimension of matrix
# tol = 1e-8				# Convergence tolerance
# mmax = n//2				# Maximum number of iterations	

# ''' Create sparse, diagonally dominant matrix A with 
#  	diagonal containing 1,2,3,...n. The eigenvalues
#     should be very close to these values. You can 
#     change the sparsity. A smaller number for sparsity
#     increases the diagonal dominance. Larger values
#     (e.g. sparsity = 1) create a dense matrix
# '''

# sparsity = 0.0001
# A = np.zeros((n,n))
# for i in range(0,n):
#     A[i,i] = i + 1 
# A = A + sparsity*np.random.randn(n,n) 
# A = (A.T + A)/2 


# k = 8					# number of initial guess vectors 
# eig = 4					# number of eignvalues to solve 
# t = np.eye(n,k)			# set of k unit vectors as guess
# V = np.zeros((n,n))		# array of zeros to hold guess vec
# I = np.eye(n)			# identity matrix same dimen as A

# # Begin block Davidson routine

# start_davidson = time.time()

# for m in range(k,mmax,k):
#     if m <= k:
#         for j in range(0,k):
#             V[:,j] = t[:,j]/np.linalg.norm(t[:,j])
#         theta_old = 1 
#     elif m > k:
#         theta_old = theta[:eig]
#     V[:,:m],R = np.linalg.qr(V[:,:m])
#     T = np.dot(V[:,:m].T,np.dot(A,V[:,:m]))
#     THETA,S = np.linalg.eig(T)
#     idx = THETA.argsort()
#     theta = THETA[idx]
#     s = S[:,idx]
#     for j in range(0,k):
#         w = np.dot((A - theta[j]*I),np.dot(V[:,:m],s[:,j])) 
#         q = w/(theta[j]-A[j,j])
#         V[:,(m+j)] = q
#     norm = np.linalg.norm(theta[:eig] - theta_old)
#     if norm < tol:
#         break

# end_davidson = time.time()

# # End of block Davidson. Print results.

# print("davidson = ", theta[:eig],";", end_davidson - start_davidson, "seconds")

# # Begin Numpy diagonalization of A

# start_numpy = time.time()

# E,Vec = np.linalg.eig(A)
# E = np.sort(E)

# end_numpy = time.time()

# # End of Numpy diagonalization. Print results.

# print("numpy = ", E[:eig],";", end_numpy - start_numpy, "seconds") 