# number of atoms : 2
# material : hydrogen
# space configuration : finite difference method (유한차분법)


import numpy as np
import random
import os

os.getcwd()

import imp
import Basis
import Density

# 2. from the veff, we can get full potential
# 3. then calculate the coefficient of electron density functions.

# imp.reload(Basis)
# imp.reload(Density)

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
K = Basis.get_K_matrix(configs)

K[0][1](r = [0.1, 0.2, 0.1])


############################
# 각 node의 인덱스 생성
############################
nx= ny =nz = 3
flattened_index_of_nodes = [None for i in range(nx * ny * nz)]
q = 0
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            flattened_index_of_nodes[q] = (i, j, k)
            q += 1


nik = Density.generate_nik_by_pw(configs)
nr = Density.cal_Nr(configs, [0.5, 0.2, 0.7], nik)


