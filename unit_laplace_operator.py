# -*- coding: utf-8 -*-

# 출처 : http://hpcschool.kr/hpcss2014/wp-content/uploads/sites/12/2014/07/FiniteDifferenceMethod.pdf

"""
    gernerate nodes
"""

# 1. solve poisson equation, get hartree potential (veff)
    # - create node

nx = ny = nz = 3

nodes = np.zeros((nx, ny, nz))

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
