'''
MAP3121 - 2022
EP3 - Modelagem de um Sistema de Resfriamento de Chips
11804608 - Arthur Pedroso Porto Belli
11965585 - Letícia Harumi Furusawa
'''
import numpy as np
from MAP3121.EP3.utils import *#contém as funções inalteradas de EPs passados

# função de excitação do sistema
def f(x, y=0):
    return 12*x*(1-x)-2

def phi_i1(y, x): #y é a variável independente, x é o array de pontos e i é o índice desejado
    return (y-x) #expressão para o trecho [x_{i-1}, x_i]

def phi_i2(y, x):
    return (x-y) # expressão de \phi para o trecho [x_i, x_{i+1}]

def funcao1(f, y, z, x): #funcao a ser integrada deduzidda na secção de exercícios do relatório
    return f(y)*phi_i1(y, x)

def funcao2(f, y, z, x): #para serem compatíveis com o código desenvolvido para o EP2, precisamos de duas vars. independentes (y, z)
    return f(y)*phi_i2(y, x)

# condutividade térmica
def k(x):
    return 1

# Soma do calor gerado pelo chip com o calor retirado pelo resfriamento, no equilibrio eles são iguais. Assim, o calor interno não se altera
def q(x):
    return 0



'''
primeira implementação, estou tentando resolver o caso mais simples
L = 1 k(x)=1 q(x)=0
'''

def fem(n):
    h = 1/(n+1) #passo dos apoios
    x = np.array([i*h for i in range(n+1)], dtype=np.double)

    #definição da matriz normal do MMQ 
    a = np.array([0 if i == 0 else (-1)*(1/h) for i in range(n)], dtype=np.double) #diagonal inferior
    b = np.array([2/h for i in range(n)], dtype=np.double) #diagonal principal
    c = np.array([0 if i == n-1 else (-1)*(1/h) for i in range(n)], dtype=np.double) #diagonal superior

    #array de termos independentes
    #putaria do krl q eu tive q fazer aqui pra funcionar, mas acho q funciona.....
    d = [(1/h)*(gauss_2((funcao1, f, 0, 1, x[i-1], x[i], 10, *x_w(10), "inf") + gauss_2(funcao2, f, 0, 1, x[i], x[i+1], 10, *x_w(10), "sup"))) for i in range(1, n)]
    d = np.array(d, dtype=np.double)

    #resolução do sistema linear
    alphas = np.array(solve_lin_sys_tridiag(*A_to_LU_tridig(a, b, c), d), dtype=np.double)


def temperatura(n): # calcula a temperatura pelo fem
    return fem(n)

def temperatura2(f,n,L): # calcula a temperatura pela equação f dada
    
    h= L/(n+1)
    x = np.array([i*h for i in range(n+1)], dtype=np.double)
    u = np.array([f(x[i]) for i in range (1,n+1)],dtype=np.double)

    return u

def erro_max(u_fem, L, n, u_teor): # devolve o erro maximo entre o valor encontrado pelo fem e o teórico

    h=L/(n+1)
    x = np.array([i*h for i in range(n+1)], dtype=np.double)
    erro = u_teor - u_fem
    erro_m = np.max(erro)

    return erro_m



def solve():
    pass

def main():
    pass


if __name__ == "__main__":
    main()

