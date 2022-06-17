# number of atoms : 2
# material : hydrogen
# space configuration : finite difference method (유한차분법)

# import packages
import numpy as np
import random




# n(~x) = ∑|ψi(~x)|^2          (2.9) , (occ에 대해 sum)

# G-space

i_max = 1; k_max = 10; G_max = 10
#######################################
# i는 occupied state
# k는 nodes의 인덱스 
# G는 reciprocal space vector, 의미가 정확히 뭐고 수학적으로 어떻게 풀어야 하지?
#######################################



C = np.random.random((i_max, k_max, G_max))


def flatten_3d_arr(nx, ny, nz): # flattened 3d array
    return np.zeros((nx - 1) * (ny - 1) * (nz - 1))

def generate_plane_wave_basis(r = None, option="functional"):

    global C, i_max, k_max, G_max
    
    Psi = [[None for k in range(k_max)] for i in range(i_max)]
    
    
    if option == "functional":
        # r 함수형
        for i in range(i_max):
            for k in range(k_max):
                def func(r, i = i, k = k):
                    global C
                    tmp = []
                    for G in range(G_max):                
                        print(C[i][k][G])
                        tmp.append(C[i][k][G] * np.exp(complex(0, (k + G) * r)))
                    return tmp
                Psi[i][k] = func
    
        return Psi
    else:
        for i in range(i_max):
            for k in range(k_max):                
                tmp = []
                for G in range(G_max):                
                    print(C[i][k][G])
                    tmp += C[i][k][G] * np.exp(complex(0, (k + G) * r))
                Psi[i][k] = tmp
        return Psi



#################################
# calculate electron density n(x)
#################################
def generate_plane_wave_basis_nik(r = None):
    
    global C, i_max, k_max, G_max
    
    Psi = [[None for k in range(k_max)]for i in range(i_max)]
    
    for i in range(i_max):
        for k in range(k_max):            
            tmp = ''
            for G in range(G_max):                
                for G_p in range(G_max):
                    # print(C[i][k][G])
                    tmp = tmp + " + " + f"{C[i][k][G]} * {C[i][k][G_p]} * np.exp(complex(0, {(k + G)} * r) - complex(0, {(k + G_p)} * r))"
            Psi[i][k] = tmp
    return Psi



nik = generate_plane_wave_basis_nik()



"""
    test nik
"""
r = 0.5
eval(nik[0][0])



# where k is node of reciperocal space
# k가 1 ~ 9이면 node 순으로는 1번이 x = y = z = 0인 위치
# reshape을 하면 첫번째 매트릭스부터 오른쪽으로 진행됨 ((0, 0), (0, 1), (0, 2), (1, 0)... )
# 이 경우 r로 변환하면?

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





"""
    gernerate nodes
"""

# 1. solve poisson equation, get hartree potential (veff)
    # - create node

n_node_of_x = n_node_of_y = n_node_of_z = 3

nodes = np.zeros((n_node_of_x, n_node_of_y, n_node_of_z))

    # - set hydrogen atom in random position
[h_x, h_y, h_z] = [random.random() * 10 for i in range(3)]

    # consider process of calculation of possion equation
    # create laplace operator, which is independent of basis functions 
    # create trial vector, we need the configuration of basis functions




# # only  5 diagonal component is non-zero for 2d
# nx = ny = nz = 3; dx = dy = dz = 1
def generate_unit_laplace_operator(nx, ny, nz, dx = 1, dy = 1, dz = 1): 
    mat = np.zeros(((nx - 1) * (ny - 1) * (nz - 1), (nx - 1) * (ny - 1) * (nz - 1)))

    # fill diagonal term
    np.fill_diagonal(mat, -4 / dx ** 2)

    # fill off-diagonal term ~ x
    rng = np.arange(mat.shape[0] - 1)
    mat[rng, rng + 1] = 1 / dx ** 2
    mat[rng, rng - 1] = 1 / dx ** 2
    
    # fill off-diagonal term ~ y
    rng = np.arange(mat.shape[0] - 1 - nx)
    mat[rng + nx, rng] = - 4 / dy ** 2
    mat[rng + nx, rng + 1] = 1 / dy ** 2
    mat[rng + nx, rng - 1] = 1 / dy ** 2
    mat[rng , rng + nx] = - 4 / dy ** 2
    mat[rng , rng + nx + 1] = 1 / dy ** 2
    mat[rng , rng + nx - 1] = 1 / dy ** 2
    
    # fill off-diagonal term ~ z
    rng = np.arange(mat.shape[0] - 1 - nx - ny)
    mat[rng + nx + ny, rng] = - 4 / dz ** 2
    mat[rng + nx + ny, rng + 1] = 1 / dz ** 2
    mat[rng + nx + ny, rng - 1] = 1 / dz ** 2
    mat[rng , rng + nx + ny] = - 4 / dz ** 2
    mat[rng , rng + nx + ny + 1] = 1 / dz ** 2
    mat[rng , rng + nx + ny - 1] = 1 / dz ** 2
    return mat

laplace_op = generate_unit_laplace_operator(nx, ny, nz)





# 2. from the veff, we can get full potential
# 3. then calculate the coefficient of electron density functions.






