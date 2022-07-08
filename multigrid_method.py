# -*- coding: utf-8 -*-

# AMG : Algebraic Multi-Grid method
# multigrid method need 2 components below
# 1. projection operator
#   construct the system on the coarse mesh from the original dense mesh, and correct the results on the dense mesh by the solution obtained from the coarse mesh

# 2. smoother
#   employed to damp out the relatively high frequency parts of the numerial error of the result, with respect to current mesh level
#   Gauss-Seidel iteration method
#   



import numpy as np
import random


# 멀티 그리드 메소드에서 phi는 just a values
"""
    A * Phi = f,
    
    여기서 Phi가 1차원 행렬이어야 함.
    f ( = nr) 도 마찬가지 
    일단 nr을 먼저 구해서 1차원으로 flatten 해보자.    
      
    
    x = [V (1, 1) V (1, 2) V (1, 3) · · · V (2, 1) V (2, 2) V (2, 3)] ^ T
    
    
    grid 구성 4 x 4
    
    
    n13- n14 - n15- n16
       
    n9 - n10 - n11- n12
        
    n5 - n6  - n7 - n8
    
    n1 - n2  - n3 - n4
    
    
    3 X 3 3차원이면?
    
    X_i = [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]


    −a/h^2(φ_(i−1,j,k)  + φ_(i+1,j,k) + φ_(i,j−1,k) + φ_(i,j+1,k) + φ_(i,j,k-1) + φ_(i,j,k+1) −  6 φi,j,k) = f_(i,j) , where i, j = 2, ..., N − 1,

    Matrix 구성은
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    .
    .
    .
    [v111, v112, v113, v121, v122, v123, v131, v132, v133, v211, v212, v213, v221, v222, v223, v231, v232, v233, v311, v312, v313, v321, v322, v323, v331, v332, v333]
    
    예를 들어 v222번 노드라고 하면 14번 인덱스에 있으며,
    
    i - 1인 v122는 4번
    i + 1은 v322는 22번
    
    j - 1: v212,  10
    j + 1: v232,  16
    
    k-1 : v221, 12
    k+1 : v223, 14
    

"""

N = 3
a = []
for i in range(N):
    for j in range(N):
        for k in range(N):
            a.append(f"v{i}{j}{k}")

x = np.array(a)

A = np.zeros((N **3, N ** 3))

eqns = []
for i in range(N):
    for j in range(N):
        for k in range(N):
            if i >= 1 and j >= 1 and k >= 1 and i <= N -1 and j <= N-1 and k <= N-1:
                eqns.append(f"-6 * v{i-1}{j}{k} + v{i-1}{j}{k} + v{i+1}{j}{k} + v{i}{j-1}{k}+ v{i}{j+1}{k} + v{i}{j}{k-1} + v{i}{j}{k+1} ")
                
                row_ind = np.where( x == f"v{i}{j}{k}")
                A[row_ind, row_ind] = -6

                ind = np.where( x == f"v{i-1}{j}{k}")
                A[row_ind, ind] = 1

                ind = np.where( x == f"v{i+1}{j}{k}")
                A[row_ind, ind] = 1

                ind = np.where( x == f"v{i}{j-1}{k}")
                A[row_ind, ind] = 1

                ind = np.where( x == f"v{i}{j+1}{k}")
                A[row_ind, ind] = 1

                ind = np.where( x == f"v{i}{j}{k-1}")
                A[row_ind, ind] = 1

                ind = np.where( x == f"v{i}{j}{k+1}")
                A[row_ind, ind] = 1















