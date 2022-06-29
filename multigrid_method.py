# -*- coding: utf-8 -*-

# AMG : Algebraic Multi-Grid method
# multigrid method need 2 components below
# 1. projection operator
#   construct the system on the coarse mesh from the original dense mesh, and correct the results on the dense mesh by the solution obtained from the coarse mesh

# 2. smoother
#   employed to damp out the relatively high frequency parts of the numerial error of the result, with respect to current mesh level
#   Gauss-Seidel iteration method
#   