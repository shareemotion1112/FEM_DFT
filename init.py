# number of atoms : 2
# material : hydrogen
# space configuration : finite difference method (유한차분법)

# import packages
import numpy as np
import random


# 1. solve poisson equation, get hartree potential (veff)
    # - create node

n_node_of_x = n_node_of_y = n_node_of_z = 3

nodes = np.zeros((n_node_of_x, n_node_of_y, n_node_of_z))

    # - set hydrogen atom in random position
[h_x, h_y, h_z] = [random.random() * 10 for i in range(3)]

    # consider process of calculation of possion equation
    # create laplace operator, which is independent of basis functions 
    # create trial vector, we need the configuration of basis functions

def generate_basis(nx, ny, nz): # flattened 3d array
    return np.zeros((nx - 1) * (ny - 1) * (nz - 1))


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



