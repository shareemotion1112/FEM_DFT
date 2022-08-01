# number of atoms : 2
# material : hydrogen
# space configuration : finite difference method (유한차분법)


import enum
import math
import numpy as np
import random
import os

os.getcwd()

import imp
import Basis
import Density
# import backup.Hatree_Potential_FDM as Hatree_Potential_FDM
import multigrid_method

def re_import_lib():
	imp.reload(Basis)
	imp.reload(Density)
	imp.reload(multigrid_method)

#################
# k-point 설정
# C, K등에 임의의 상수 배정 
#################
i_max = G_max = 1 # index i는 occupied state, G는 고체내 주기성을 표현하는 reciprocal vector
kx = ky = kz = 9 # k는 k-points
size_of_mesh = 10 # 나중에 메인 설정에서 받아와야 함

Div = 100

k_m = [None for i in range(kx * ky * kz)]

for i in range(len(k_m)):
    k_m[i] = [np.random.random() * size_of_mesh, np.random.random() * size_of_mesh, np.random.random() * size_of_mesh / Div ] 

G_m = [[0, 0, 0]]

k_max = len(k_m)
C = np.random.random((i_max, k_max, G_max)) / Div
configs = {'i_max' : i_max, 'k_max' : len(k_m), 'G_max' : G_max, 'C' : C, 'k_m' : k_m, 'G_m' : G_m}

# base = Basis.generate_plane_wave_basis(configs)

############################
# 각 node의 인덱스 생성
# n_{ik} 계산
# laplace operator 생성
############################

nx = ny = nz = 6; 
h = 1e-10 # 옹스트롱 단위로 진행.

flattened_index_of_nodes = [None for i in range(nx * ny * nz)]
q = 0
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            flattened_index_of_nodes[q] = (i, j, k)
            q += 1

nik = Density.generate_nik_by_pw(configs)


# 라플라서 operator 만들때도 경게조건을 넣을 수 있지 않나??  ------------ 2022.08.01
#  출처 : 42번식  in   https://my.ece.utah.edu/~ece6340/LECTURES/Feb1/Nagel%202012%20-%20Solving%20the%20Generalized%20Poisson%20Equation%20using%20FDM.pdf 
A, ind_vector = multigrid_method.create_laplace_matrix(nx, ny, nz)

f = np.zeros(len(ind_vector))
for i, ind in enumerate(ind_vector):
    x = int(ind[1]) * h
    y = int(ind[2]) * h
    z = int(ind[3]) * h
    r = [x, y, z]    
    f[i] = Density.cal_Nr(configs, r, nik)

HF_potential = np.inner(np.linalg.inv(A), -1 * f) 
print(f"HF_potential.shape : {HF_potential.shape}, ind_vector.shape : {ind_vector.shape}")

############################################################################
# HF potential 계산한거 검증 필요 ----------------------------------  2022.07.08
# 경계 조건을 어딘가에 넣어야 할 것 같은데?? ----------------------------- 2022.08.01
# 경계조건은 핵에 의한 coulomb iteraction에 의한 효과를 중첩하면 됨.
# HF potential에 핵과의 전자기력에 의한 coulomb interaction 포텐셜을 더하고 
# exchange correlation에 의한 포텐셜까지 합산해서 v_{eff}를 구하자.
############################################################################


# 쿨롱 interaction 계산
# 수소 원자 위치 설정
# grid 구조 분석
print(f"nx, ny, nz, h : {nx} {ny} {nz} {h}")

# 수소 원자 갯수 설정
N_Hatoms = 2
charges = [1, 1]


# 수소원자 position 설정
init_poss = []
for i in range(N_Hatoms):
	init_pos = [np.random.random() * nx * h for i in range(3)]
	init_poss.append(init_pos)

electron_permittivity = 8.854 * 1e-12   # [F/m]
dielectric_constant = 1

v_columb = np.zeros(len(ind_vector))
for pos in init_poss:
	for i, ind in enumerate(ind_vector):
		x = float(ind[1]) * h
		y = float(ind[2]) * h
		z = float(ind[3]) * h
		r = np.sqrt(np.sum([(x - pos[0])**2, (y - pos[1])**2, (z - pos[2])**2]))
		print(r)
		v_columb[i] = 1 / (4 * math.pi * electron_permittivity * r)



# exchange correlation 추가 예정...



# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
import matplotlib.pyplot as plt    

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

 



# 2. from the veff, we can get full potential
# 3. then calculate the coefficient of electron density functions.

