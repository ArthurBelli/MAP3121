'''
MAP3121 - 2022
EP3 - Modelagem de um Sistema de Resfriamento de Chips
11804608 - Arthur Pedroso Porto Belli
11965585 - Letícia Harumi Furusawa
'''
import numpy as np
from MAP3121.EP3.utils import *#contém as funções inalteradas de EPs passados

# função de excitação do sistema
def f(x, y):
    return 12*x*(1-x)-2

# condutividade térmica
def k(x):
    return 1

# Soma do calor gerado pelo chip com o calor retirado pelo resfriamento, no equilibrio eles são iguais. Assim, o calor interno não se altera
def q(x):
    return 0





def fem(L, n):
    h = L/(n+1)
    x = np.array([i*h for i in range(n+1)], dtype=np.double)

    #definição da matriz normal do MMQ 
    a = np.array([0 if i == 0 else (-1)*(1/h) for i in range(n)], dtype=np.double) #diagonal inferior
    b = np.array([2/h for i in range(n)], dtype=np.double) #diagonal principal
    c = np.array([0 if i == n-1 else (-1)*(1/h) for i in range(n)], dtype=np.double) #diagonal superior

    # vetor de termos independentes
    #imagino que tenha q usar expressão lambda, mas n sei se a implementação ta certa, vou dar uma pesquisada no fim de semana
    d = [(1/h)*(gauss_2((lambda x:f(x))*((lambda x:x)-x[i-1]), 0, 1, x[i-1], x[i], 10, *x_w(10))+gauss_2((lambda x: f(x))*(x[i+1]- (lambda x:x), 0, 1, x[i], x[i+1], 10, *x_w(10)))) for i in range(1, n)]


def solve():
    pass

def main():
    pass


if __name__ == "__main__":
    main()

