# -*- coding: utf-8 -*-
# 출처 : https://www.notion.so/basis-set-72f97e7e902e48f2866d63a4f39b1a4e

# 출처2 : 갤러킨 방법, 가중잔여법(최소제곱법, 갤러킨방법)
#   https://m.blog.naver.com/PostView.naver?isHttpsRedirect=true&blogId=mykepzzang&logNo=221113579060
#   https://m.blog.naver.com/PostView.naver?blogId=mykepzzang&logNo=221114379052&targetKeyword=&targetRecommendationCode=1
#       라플라스 방정식의 해를 여러 함수의 중첩인 근사해로 가정하고
#           d^2u/dx^2 + p(x) = 0. 여기 u는 exact solution
#           d^2u'/dx^2 + p(x) = R(x),   여기 u'는 approximate solution
#               u` = simga c_i * phi_i
#           integral R * phi_i * dx = integral {d^2u'/dx^2 + p(x)} * phi_i * dx = 0
#           
#           부분적분을 진행하면,
# 
#           [K][c] = [F], 
#           [c] = [K]^(-1) [F]
#
#           [K] = integral(0 ~ 1) dphi_i / dx * dphi_j / dx * dx
#           [F] = integral(0 ~ 1) p(x) * [phi(x)] * dx + du`/dx(1) [phi(1)] * du'/dx(0) [phi(0)], 여기서 0과 1은 geometry(mesh 및 노드)의 최소, 최대값

# import packages
import numpy as np
import random




# n(~x) = ∑|ψi(~x)|^2          (2.9) , (occ에 대해 sum)

# G-space

i_max = 1; k_max = 10; G_max = 10
#######################################
# i는 occupied state
# k는 nodes의 인덱스, discrete fourier transform 처럼 k가 각 축의 node 개수가 되는 듯 
# G는 reciprocal space vector, 의미가 정확히 뭐고 수학적으로 어떻게 풀어야 하지?
#######################################


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



##################
# C 초기화
##################
C = np.random.random((i_max, k_max, G_max))


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

basis = generate_plane_wave_basis()

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

