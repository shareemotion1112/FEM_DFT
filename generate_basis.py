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


#######################################
# 출처 : https://www.archer.ac.uk/training/course-material/2014/04/PMMP_UCL/Slides/castep_1.pdf
# i는 occupied state
# k는 nodes의 인덱스, discrete fourier transform 처럼 k가 각 축의 node 개수가 되는 듯 
#   k-point 설정을 통해 얻어야 하는 값임. 아래 설명(k-point 설)을 참조하여 생성해야 함.
# G는 reciprocal space vector, 의미가 정확히 뭐고 수학적으로 어떻게 풀어야 하지?
# G는 sysmetric한 방향으로의 vector를 의미하는 듯, 비정질 구조에서는 일단 제외 가능해보임.
#   일단, fourier 변환 후 발생하는 frequency의 일종이므로 
#        Psi_(i, k) (r) = sum_(G) C_(i,k) * exp(i(k + G) * r), where Psi는 위 설명의 phi이자 trial functions
#   foureir 변환을 시행할 수 밖에 없는 듯??

# k point 설정 : https://www.materialssquare.com/blog/2-convergence-test-k-points-optimization-for-silicon-bulk-2-ko
# Periodic boundary cell에서 파동 함수는 Bloch theorem을 이용하여 나타낼 수 있습니다.[1]
# 이때 파동함수는 Brillouin zone에서의 벡터인 k 에 대해 서술되는데, 
# Brillouin zone이란 역격자 공간에서 수직이등분선으로 이루어진 가장 작은 다면체를 말하며, 
# 원래 격자와 같은 point symmetry를 가집니다.
# Brillouin zone에 존재하는 수많은 wave-vector k 모두에 대해 해를 구하면 시스템을 가장 잘 서술할 수 있지만, 
# 현실적으로 어려운 일이므로 k-point sampling을 통해 그 중 일부를 골라내어 계산을 수행해야 합니다. 
# 통계에서 오차를 줄이기 위해 표본 집단을 여러 방식으로 선택하는 것처럼, 
# k-point를 sampling하는 방법에는 여러 가지가 있습니다. 
# 그 중에서 가장 일반적이며 널리 사용되는 방법이 Monkhorst-Pack sampling입니다.
# 이 방법은 전 공간에 대해 symmetry를 고려하여 균등한 간격으로 sampling하는 방법으로, 
# 특별한 목적이 없다면 이 옵션을 사용하는 것이 적절합니다. 
# Quantum espresso에서는 automatic 옵션을 선택하면 Monkhorst-Pack 방식으로 k-point sampling을 할 수 있습니다.
#######################################



#################
# k-point 설정
#################
kx = ky = kz = 9
size_of_mesh = 10 # 나중에 메인 설정에서 받아와야 함

k_m = [None for i in range(kx * ky * kz)]

for i in range(len(k_m)):
    k_m[i] = [np.random.random() * size_of_mesh, np.random.random() * size_of_mesh, np.random.random() * size_of_mesh]

G_m = [[0, 0, 0]]
i_max = 1; k_max = len(k_m); G_max = 1


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
                    for g in range(G_max):                
                        print(C[i][k][G])
                        tmp.append(C[i][k][G] * np.exp(complex(0, (k_m[k] + G_m[g]) * r)))
                    return tmp
                Psi[i][k] = func
    
        return Psi
    else:
        for i in range(i_max):
            for k in range(k_max):                
                tmp = []
                for g in range(G_max):                
                    print(C[i][k][G])
                    tmp += C[i][k][G] * np.exp(complex(0, (k_m[k] + G_m[g]) * r))
                Psi[i][k] = tmp
        return Psi

basis = generate_plane_wave_basis()

#################################
# calculate electron density n(x)
#################################
def generate_plane_wave_basis_nik(r = None):
    
    global C, i_max, k_max, G_max, k_m, G_m
    
    Psi = [[None for k in range(k_max)]for i in range(i_max)]
    
    for i in range(i_max):
        for k in range(k_max):            
            tmp = ''
            for g in range(G_max):                
                for g_p in range(G_max):
                    print(C[i][k][g])
                    tmp = tmp + " + " + f"{C[i][k][g]} * {C[i][k][g_p]} * \
                    np.exp(complex(0, { k_m[k][0] + G_m[g][0] } * r[0] ) - complex(0, { k_m[k][0] + G_m[g_p][0] } * r[0] ) \
                            + complex(0, { k_m[k][1] + G_m[g][1] } * r[1] ) - complex(0, { k_m[k][1] + G_m[g_p][1] } * r[1] ) \
                            + complex(0, { k_m[k][2] + G_m[g][1] } * r[2] ) - complex(0, { k_m[k][2] + G_m[g_p][2] } * r[2] ))"
            Psi[i][k] = tmp
    return Psi

a = f"{C[i][k][g]} * {C[i][k][g_p]} * np.exp(complex(0, { k_m[k][0] + G_m[g][0] } * {r[0]} ) - complex(0, { k_m[k][0] + G_m[g_p][0] } * {r[0]} ) \
        + complex(0, { k_m[k][1] + G_m[g][1] } * {r[1]} ) - complex(0, { k_m[k][1] + G_m[g_p][1] } * {r[1]} ) \
        + complex(0, { k_m[k][2] + G_m[g][1] } * {r[2]} ) - complex(0, { k_m[k][2] + G_m[g_p][2] } * {r[2]} ) ) "
eval(a)

nik = generate_plane_wave_basis_nik()

###################################
# n(r) = 1/N_k * sum_(k, i)(n_ik)
###################################

"""
    test nik
"""
r = [0.5, 0.4, 0.2]
eval(nik[0][0])



# where k is node of reciperocal space
# k가 1 ~ 9이면 node 순으로는 1번이 x = y = z = 0인 위치
# reshape을 하면 첫번째 매트릭스부터 오른쪽으로 진행됨 ((0, 0), (0, 1), (0, 2), (1, 0)... )
# 이 경우 r로 변환하면?

