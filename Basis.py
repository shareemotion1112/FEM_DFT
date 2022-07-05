# -*- coding: utf-8 -*-
# 출처 : https://www.notion.so/basis-set-72f97e7e902e48f2866d63a4f39b1a4e



# import packages
import numpy as np

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


