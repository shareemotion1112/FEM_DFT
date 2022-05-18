# number of atoms : 2
# material : hydrogen
# space configuration : finite difference method

# import packages
import numpy as np
import random


# 1. solve poisson equation, get hartree potential (veff)
    # - create node

n_node_of_x = 10
n_node_of_y = 10
n_node_of_z = 10

nodes = np.zeros((n_node_of_x, n_node_of_y, n_node_of_z))

    # - set hydrogen atom in random position
[h_x, h_y, h_z] = [random.random() * 10 for i in range(3)]

    # consider process of calculation of possion equation
    # create laplace operator, which is independent of basis functions 
    # create trial vector, we need the configuration of basis functions

def generate_basis(nx, ny, nz): # flattened 3d array
    return np.zeros((nx - 1) * (ny - 1) * (nz - 1))


# nx = ny = nz = 10; dx = dy = dz = 1
def unit_laplace_operator(nx, ny, nz, dx, dy, dz): 
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
    
    # fill off-diagonal term ~ y
    rng = np.arange(mat.shape[0] - 1 - nx - ny)
    mat[rng + nx + ny, rng] = - 4 / dz ** 2
    mat[rng + nx + ny, rng + 1] = 1 / dz ** 2
    mat[rng + nx + ny, rng - 1] = 1 / dz ** 2
    mat[rng , rng + nx + ny] = - 4 / dz ** 2
    mat[rng , rng + nx + ny + 1] = 1 / dz ** 2
    mat[rng , rng + nx + ny - 1] = 1 / dz ** 2
    return mat


# # only  5 diagonal component is non-zero for 2d
# def laplace_operator(u_matrix, dx, dy, dz):
#     result = np.zeros_like(u_matrix)

#     for i in range(1, n_node_of_x - 1):
#         for j in range(1, n_node_of_y - 1):
#             for k in range(1, n_node_of_z - 1):
#                 f_x = (result[i + 1, j, k] + result[i - 1, j, k] - 2 * result[i, j, k]) / dx ** 2
#                 f_y = (result[i, j + 1, k] + result[i, j - 1, k] - 2 * result[i, j, k]) / dy ** 2
#                 f_z = (result[i, j, k + 1] + result[i, j, k - 1] - 2 * result[i, j, k]) / dz ** 2
#                 result[i, j, k] = f_x + f_y + f_z
#     return result

# 2. from the veff, we can get full potential
# 3. then calculate the coefficient of electron density functions.



