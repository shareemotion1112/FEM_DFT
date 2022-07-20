# number of atoms : 2
# material : hydrogen
# space configuration : finite difference method (유한차분법)


import enum
from matplotlib import projections
import numpy as np
import random
import os

os.getcwd()

import imp
import Basis
import Density
import backup.Hatree_Potential_FDM as Hatree_Potential_FDM
import multigrid_method

# 2. from the veff, we can get full potential
# 3. then calculate the coefficient of electron density functions.

# imp.reload(Basis)
# imp.reload(Density)
# imp.reload(multigrid_method)

#################
# k-point 설정
#################
i_max = G_max = 1
kx = ky = kz = 9
size_of_mesh = 10 # 나중에 메인 설정에서 받아와야 함

k_m = [None for i in range(kx * ky * kz)]

for i in range(len(k_m)):
    k_m[i] = [np.random.random() * size_of_mesh, np.random.random() * size_of_mesh, np.random.random() * size_of_mesh]

G_m = [[0, 0, 0]]

k_max = len(k_m)
C = np.random.random((i_max, k_max, G_max))
configs = {'i_max' : i_max, 'k_max' : len(k_m), 'G_max' : G_max, 'C' : C, 'k_m' : k_m, 'G_m' : G_m}


base = Basis.generate_plane_wave_basis(configs)
# K = Hatree_Potential_FDM.get_K_matrix(configs) # 적분하기 전임.. 적분해야함
# F1 = Hatree_Potential_FDM.get_F_integral_part_matrix(configs, r=[0.1, 0.1, 0,1])


############################
# 각 node의 인덱스 생성
############################
nx = ny = nz = 6; h = 1
flattened_index_of_nodes = [None for i in range(nx * ny * nz)]
q = 0
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            flattened_index_of_nodes[q] = (i, j, k)
            q += 1


nik = Density.generate_nik_by_pw(configs)
nr = Density.cal_Nr(configs, [0.5, 0.2, 0.7], nik)

A, ind_vector = multigrid_method.create_laplace_matrix(nx, ny, nz)

f = np.zeros(len(ind_vector))
for i, ind in enumerate(ind_vector):
    x = int(ind[1]) * h
    y = int(ind[2]) * h
    z = int(ind[3]) * h
    r = [x, y, z]    
    f[i] = Density.cal_Nr(configs, r, nik)

HF_potential = np.inner(np.linalg.inv(A), -1 * f) 

# HF potential 계산한거 검증 필요 2022.07.08
HF_potential.shape
ind_vector.shape
    
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt    
from matplotlib import cm


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

size = ind_vector.shape[0]

x = np.zeros(size)
y = np.zeros(size)
z = np.zeros(size)
c = np.zeros(size)

for i, ind in enumerate(ind_vector):
    x[i] = int(ind[1])
    y[i] = int(ind[2])
    z[i] = int(ind[3])
    c[i] = HF_potential[i]


c_norm = (c - np.min(c)) / (np.max(c) - np.min(c))
c_norm = [round(c, 2) for c in c_norm]
mask = c_norm > 0.1

img  = ax.scatter(x, y, z, c = c_norm, s = 20)
plt.tight_layout()
fig.colorbar(img)
plt.show()

 