'''
MAP3121 - 2022
EP3 - Modelagem de um Sistema de Resfriamento de Chips
11804608 - Arthur Pedroso Porto Belli
11965585 - Letícia Harumi Furusawa
'''

import numpy as np
import matplotlib.pyplot as plt
from utils import * #contém as funções inalteradas de EPs passados


def zero(x):
    return 0

# condutividade térmica

def k(x):
    return 3.6

def que_const(x, Qp0=400, Qm0=200):
    return Qp0-Qm0

def que(x, Qp0=600, Qm0=300, sigma=0.1, theta=1, L=1):
    q = Qp0*np.exp(-(x-L/2)**2/sigma**2) - Qm0*(np.exp(-(x)**2/theta**2) + np.exp(-((x-L)**2/theta**2)))
    return q

# função de excitação do sistema
def excitacao(x):
    return 12*x*(1-x)-2

def phi(x, apoios, h, i): #funções da base do espaço gerado para o método dos elementos finitos
    if apoios[i-1] <= x <= apoios[i]: #primeira metade do spline é crescente
        return (x-apoios[i-1])/h
    elif apoios[i] <= x <= apoios[i+1]: #segunda metade do spline é decrescente
        return (apoios[i+1]-x)/h
    return 0 # fora do intervalo

def debug_phi(apoios, h, n):
    for k in range(1, n+1):
        print(f"phi_{k}_subida = (x-{apoios[k-1]}/{h}), para x in [{apoios[k-1]}, {apoios[k]}]")
        print(f"phi_{k}_descida = ({apoios[k+1]} - x/{h}), para x in [{apoios[k]}, {apoios[k+1]}]")

def calc_v(alphas, x, apoios, h, n):
    temp = float(sum([alphas[i-1]*phi(x, apoios, h, i) for i in range(1, n+1)]))
    return temp

def fem(n, L, k, f):
    '''
    retorna a solução do problema normalizado em x in [0,1]
    n -> número de apoios
    L -> comprimento do chip
    k -> condutividade térmica do material
    f -> função de excitação do sistema
    '''
    h = L/(n+1) #passo dos apoios

    #print(f"h = {h}")

    apoios = [i*h for i in range(n+2)]
    #print("apoios não normalizados: ")
    #print(apoios)

    apoios_barr = [apoios[i]/L for i in range(n+2)]
    #print("apoios normalizados: ")
    #print(apoios_barr)
    #print("\n")

    #debug_phi(apoios_barr, h, n)
    #definição da matriz normal do MMQ 

    #diagonal inferior
    #print("Diagonal inferior")
    a = [0 if i == 0 else (-1)*(1/(L*(h**2)))*gauss(k, apoios_barr[i], apoios_barr[i+1], 2, *x_w(2)) for i in range(n)] 
    #print(a)
    a = np.array(a, dtype=np.double)

    #diagonal principal
    #print("Diagonal principal")
    b = [(1/(L*(h**2)))*(gauss(k, apoios_barr[i-1], apoios_barr[i], 2, *x_w(2)) + gauss(k, apoios_barr[i], apoios_barr[i+1], 2, *x_w(2))) for i in range(1, n+1)] 
    #print(b)
    b = np.array(b, dtype=np.double)

    #diagonal superior 
    #print("Diagonal superior")
    c = [0 if i == n-1 else (-1)*(1/(L*(h**2)))*gauss(k, apoios_barr[i], apoios_barr[i+1], 2, *x_w(2)) for i in range(n)]
    #print(c)
    c = np.array(c, dtype=np.double)

    #array de termos independentes
    #print("Termos independentes")
    d = [(gauss((lambda nos: f(nos)*phi(nos, apoios, h, i)), apoios[i-1], apoios[i], 2, *x_w(2)) + gauss((lambda nos: f(nos)*phi(nos, apoios, h, i)), apoios[i], apoios[i+1], 2, *x_w(2))) for i in range(1, n+1)]
    #print(d)
    
    #resolução do sistema linear
    #print("alphas")
    alphas = solve_lin_sys_tridiag(*A_to_LU_tridig(a, b, c), d)
    #print(alphas)

    eixox = np.linspace(0, L, 100)
    v = [0]*len(eixox)
    for k in range(len(eixox)):
        v[k] = calc_v(alphas, eixox[k], apoios, h, n)
    return v #vetor solução normalizado

def solve_generic(n, L, k, f, a, b): # calcula a temperatura pela equação f dada
    '''
    n -> número de apoios a serem utilizados
    L -> comprimento do chip
    k -> condutividade térmica do materia
    f -> função de excitação térmica do chip
    a -> condição de contorno no extremo 0
    b -> condição de contorno no extremo L
    '''

    eixox = np.linspace(0, L, 100)
    v = fem(n, L, k, f)
    u = [0]*len(eixox)
    for i in range(len(eixox)):
        u[i] = (v[i] + a + (b-a)*eixox[i]/L)
    return u #vetor solução individualizado

def erro_max(u_teor, u_fem): # devolve o erro maximo entre o valor encontrado pelo fem e o teórico
    erro = u_teor - u_fem
    erro_m = np.max(erro)
    return erro_m

'''
parte que está dando problema
'''
def k_mix(x, L=1, d=0.05):
    if 0 <= x <= L/2 - d:
        return 60
    elif L/2 - d <= x <= L/2 + d:
        return 3.6
    elif L/2 + d <= x <= L:
        return 60
    return 0


def mat_varia(q,L,n):
  
  material21 = solve_generic(n, L, k_mix, q, 0, 0)
  return material21
'''
fim da parte que está dando problema
'''
def main():
  n = 50
  L = 1
  #u = solve_generic(n, L, k, que_const, 0, 0)
  #print(u)
  u = mat_varia(que_const, L, n)
  print(u)

    
if __name__ == "__main__":
    main()


