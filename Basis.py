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


##################
# C 초기화
##################


def generate_plane_wave_basis(config, r = None):

    i_max = config['i_max']
    k_max = config['k_max']
    G_max = config['G_max']
    k_m = config['k_m']
    G_m = config['G_m']
    C = config['C']
    
    
    Psi = [[None for k in range(k_max)] for i in range(i_max)]
    
    for i in range(i_max):
        for k in range(k_max):
            def func(r, i = i, k = k, C = C):
                tmp = []
                for g in range(G_max):                
                    # print({C[i][k][g]})
                    tmp.append({C[i][k][g]} * np.exp(complex(0, (np.array({k_m[k]}) + np.array({G_m[g]})) * np.array(r))))
                return tmp
            Psi[i][k] = func

    return Psi


# 1차원 수식
#   phi_i의 i인덱스는 trial function의 index Psi = sum ( Psi_0 + Psi_1 + ... Psi_N) 에서 N, i_max(occupied number)로 봐도 되는가? 
#   K = integral_(0, 1) * dphi_i / dx * dphi_j/ dx * dx
#   K = integral_(0, 1) grad(phi_i) * grad(phi_j) * dr_vec
#   K = integral_x(0, 1)integral_y(0, 1)integral_z(0, 1) * {dphi_i / dx * x_hat  + dphi_i / dy * y_hat + dphi_i / dz * z_hat} * {dphi_j / dx * x_hat  + dphi_j / dy  * y_hat + dphi_j / dz * z_hat} * d(x_hat + y_hat + z_hat)


def get_K_matrix(config, r=None):
    i_max = config['i_max']
    k_max = config['k_max']
    G_max = config['G_max']
    k_m = config['k_m']
    G_m = config['G_m']
    C = config['C']
    
    K = [[None for k in range(k_max)] for i in range(k_max)]

    for k in range(k_max):
        for k_p in range(k_max):
            # print(k_m[k_p])
            def func(r, k = k, k_p = k_p, C = C):
                tmp = None
                for i in range(i_max):
                    for g in range(G_max):                
                        # print(C[i][k][g])
                        if tmp == None:
                            # 미분해서 i(kx + ky + kz)가 앞에 곱해짐..
                            tmp = complex(0, np.sum(np.array(k_m[k]) + np.array(G_m[g]))) * C[i][k][g] * np.exp(  complex( 0, np.inner( (np.array(k_m[k]) + np.array(G_m[g])), np.array(r) ) )  )  * \
                            complex(0, np.sum(np.array(k_m[k_p]) + np.array(G_m[g]))) * C[i][k_p][g] * np.exp(  complex( 0, np.inner( (np.array(k_m[k_p]) + np.array(G_m[g])), np.array(r) ) )  )                        
                        else:
                            tmp += complex(0, np.sum(np.array(k_m[k]) + np.array(G_m[g]))) * C[i][k][g] * np.exp(  complex( 0, np.inner( (np.array(k_m[k]) + np.array(G_m[g])), np.array(r) ) )  )  * \
                            complex(0, np.sum(np.array(k_m[k_p]) + np.array(G_m[g]))) * C[i][k_p][g] * np.exp(  complex( 0, np.inner( (np.array(k_m[k_p]) + np.array(G_m[g])), np.array(r) ) )  )      
                return tmp
            K[k][k_p] = func

    return K









