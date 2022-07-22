'''
Como dependemos de funções implementadas nos outros EPs, 
as funções necessárias foram copiadas e, a depender do caso, adaptadas para esse arquivo a 
fim de facilitar a importação das funcionalidades sem a necessidade 
de se ter os arquivos originais no mesmo diretório
'''

import numpy as np
from math import *

# Funções do EP1
def A_to_LU_tridig(a, b, c): #recebe as diagonais
    dim = b.shape[0] #número de elementos na diagonal principal é a dimmensão da matriz
    l_ip1 = np.zeros(dim, dtype=np.double) #l_i+1
    u_ii = np.zeros(dim, dtype=np.double) #u_ii
    u_ip1 = c.copy() #u_i+1

    u_ii[0] = b[0]
    for i in range(1, dim):
        l_ip1[i] = a[i]/u_ii[i-1] #multiplicadores
        u_ii[i] = b[i]-l_ip1[i]*c[i-1]

    return u_ii, u_ip1, l_ip1

#recebe as matrizes L, U já armazenadas de forma otimizada e retorna os valores que resolvem Ax = d
def solve_lin_sys_tridiag(u, u_ip1, l, d):
    #Ly = d
    dim = u.shape[0]
    y = np.zeros((dim, 1), dtype=np.double) #vetor coluna
    y[0][0] = d[0]

    for i in range(1, dim):
        y[i][0] = d[i] - l[i]*y[i-1][0]
    
    #Ux = y
    x = np.zeros((dim, 1), dtype=np.double) #vetor coluna
    x[dim-1] = y[dim-1][0]/u[dim-1]

    for i in range(dim-1, -1, -1):
        x[i-1] = (y[i-1][0]-u_ip1[i-1]*x[i])/u[i-1]

    return x

## Funções do arquivo EP2.py

def x_w(n):
    if n == 2:
        x = [-0.5773502691896257645091488, 0.5773502691896257645091488]
        w = [1.0, 1.0]
    return x, w

# formula da gauss com mudança de variável
def find_lin_scale_transp(a, b):
    '''
    função retorna os fatores de correção de intervalos [a,b] para [-1,1]
    '''
    return np.double((b-a)/2), np.double((a+b)/2)

def gauss(f, lim_inf, lim_sup, n, nos, pesos):
    '''
    integração simples no intervalo [a,b] da função f(x) com n nós (x) e pesos (w)
    '''
    linear_scaling, linear_transposition = find_lin_scale_transp(lim_inf, lim_sup)
    sum = 0
    v = [linear_scaling*pesos[i] for i in range(n)] # pesos adaptados
    y = [(linear_scaling*nos[i]+linear_transposition) for i in range(n)] # nós adaptados
    for i in range (n): # soma iterada da integral
        sum += v[i]*f(y[i])

    return sum


def gabarito(x):
    return (x**2)*(1-x)**2