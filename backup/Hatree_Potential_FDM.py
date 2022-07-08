# 출처 : 갤러킨 방법, 가중잔여법(최소제곱법, 갤러킨방법)
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



# 1차원 수식
#   phi_i의 i인덱스는 trial function의 index Psi = sum ( Psi_0 + Psi_1 + ... Psi_N) 에서 N, i_max(occupied number)로 봐도 되는가? 
#   K = integral_(0, 1) * dphi_i / dx * dphi_j/ dx * dx
#   K = integral_(0, 1) grad(phi_i) * grad(phi_j) * dr_vec
#   K = integral_x(0, 1)integral_y(0, 1)integral_z(0, 1) * {dphi_i / dx * x_hat  + dphi_i / dy * y_hat + dphi_i / dz * z_hat} * {dphi_j / dx * x_hat  + dphi_j / dy  * y_hat + dphi_j / dz * z_hat} * d(x_hat + y_hat + z_hat)

import numpy as np
import Density
#######################
# create trial function

# u' = sum_i^N c_i * np.exp(complex(0, k_i * r))
# 너무 복잡하니까 1차 함수의 결합으로 계산??
# u'(r) = {c_1x*x + c_2x*x^2} * x_hat + {c_1y*y + c_2y*y^2} * y_hat + { c_1z * z + c_2z * z^2 } * z_hat
# 일단 x, y, z 축 각각 구해보자.


#######################

number_of_trial_functions = 5

D = [Dx,  Dy,  Dz] = [10, 10, 10] # 총길이
dr = 0.1 # 적분 구간




def get_K_matrix(r=None):    
    global D, dr
    K = [[[None for j in range(1, number_of_trial_functions)] for k in range(1, number_of_trial_functions)] for i in range(3)]    # =>  [Kx, Ky, Kz]
    
    for i in range(1, number_of_trial_functions):    
        for j in range(1, number_of_trial_functions):    
            for k in range(3):                
                tmp = 0
                for dd in np.arange(0, D[k], dr):
                    tmp += i * dd ** (i-1) * j * dd * (j-1) * dr
                K[k][i-1][j-1] = tmp
    return K

# test
# K = get_K_matrix()
# K[0]

# 경계 조건 필요
# u(0) = 0, du/dr(D) = 1
# F도 마찬가지로 [Fx, Fy, Fz] 로 구분 --------------------------> 문제점이 여기서 발생... 경계 조건이 있어야 함...
# 무엇보다 3차원 확장이 이상하네.. [Fx, Fy, Fz] 이렇게 푸는게 아닌거 같은데...

def get_F_matrix(config, nik):    
    global D, dr
    F = [[None for k in range(1, number_of_trial_functions)] for i in range(3)]    # =>  [Kx, Ky, Kz]
    
    for i in range(1, number_of_trial_functions):        
        for k in range(3):                
            tmp = 0
            r = [0, 0, 0]
            for dd in np.arange(0, D[k], dr):
                r[k] = dd
                tmp += Density.cal_Nr(config, r, nik) * dd ** (i)* dr
            tmp += D ** {i} - 0

            F[k][i] = tmp
    return F









# 1치원 trial function을 expression으로 생성?
# phi_i = ''
# for i in range(1, number_of_trial_functions):    
#     if i < (number_of_trial_functions - 1):
#         phi_i += f"x**{i}"  + ' + ' 
#     else:
#         phi_i += f"x**{i}"

# dphi_i = ''
# for i in range(1, number_of_trial_functions):    
#     if i < (number_of_trial_functions - 1):
#         dphi_i += f"{i} * x**{i-1}"  + ' + ' 
#     else:
#         dphi_i += f"{i} * x**{i-1}"

# def get_F_integral_part_matrix(config, r=None):
#     i_max = config['i_max']
#     k_max = config['k_max']
#     G_max = config['G_max']
#     k_m = config['k_m']
#     G_m = config['G_m']
#     C = config['C']

#     nik = Density.generate_nik_by_pw(config)

#     def nr(r, config, nik):
#         return Density.cal_Nr(config, r, nik)


#     F1 = [None for i in range(k_max)]

#     for k in range(k_max):
#             # print(k_m[k_p])
#         def func(r, k = k, C = C):
#             tmp = None
#             for i in range(i_max):
#                 for g in range(G_max):                
#                     # print(C[i][k][g])
#                     if tmp == None:
#                         # 미분해서 i(kx + ky + kz)가 앞에 곱해짐..
#                         tmp = nr(r, config, nik) * C[i][k][g] * np.exp(  complex( 0, np.inner( (np.array(k_m[k]) + np.array(G_m[g])), np.array(r) ) )  )                        
#                     else:
#                         tmp += nr(r, config, nik) * C[i][k][g] * np.exp(  complex( 0, np.inner( (np.array(k_m[k]) + np.array(G_m[g])), np.array(r) ) )  )     
#             return tmp
#         F1[k] = func
#     return F1



# def get_F_second_part_matrix(config, r=None):
#     i_max = config['i_max']
#     k_max = config['k_max']
#     G_max = config['G_max']
#     k_m = config['k_m']
#     G_m = config['G_m']
#     C = config['C']

#     F2 = [None for i in range(k_max)]

#     for k in range(k_max):
#         def func(r, k = k, C = C):
#             tmp = None
#             for i in range(i_max):
#                 for g in range(G_max):                
#                     # print(C[i][k][g])
#                     if tmp == None:
#                         # 미분해서 i(kx + ky + kz)가 앞에 곱해짐..
#                         tmp = complex(0, np.sum(np.array(k_m[k]) + np.array(G_m[g]))) * C[i][k][g] * np.exp(  complex( 0, np.inner( (np.array(k_m[k]) + np.array(G_m[g])), np.array(r) ) )  )                        
#                     else:
#                         tmp += complex(0, np.sum(np.array(k_m[k]) + np.array(G_m[g]))) * C[i][k][g] * np.exp(  complex( 0, np.inner( (np.array(k_m[k]) + np.array(G_m[g])), np.array(r) ) )  )     
#             return tmp
#         F2[k] = func

#     return F2