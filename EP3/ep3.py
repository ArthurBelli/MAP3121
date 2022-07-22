'''
MAP3121 - 2022
EP3 - Modelagem de um Sistema de Resfriamento de Chips
11804608 - Arthur Pedroso Porto Belli
11965585 - Letícia Harumi Furusawa
'''

import numpy as np
import matplotlib.pyplot as plt
from utils import * #contém as funções inalteradas de EPs passados


# condutividade térmica

def k(x):
    return 1

def k36(x):
    return 3.6

def k60(x):
    return 60

'''
funções do teste adicional
'''
def kex(x):
    return np.exp(x)

def qex(x):
    return (np.exp(x)+1)

def gabex():
    u = []
    eixo = np.linspace(0, 1, 100)
    for x in eixo:
        u.append((x-1)*(np.exp(-x)-1))
    return u

'''
fim do teste adicional
'''


def k_mix(x, L=1, d=0.05): #k(x) que resolve exercício 4.4
    if (0 <= x <= L/2 - d):
        return 60
    elif (L/2 - d <= x <= L/2 + d):
        return 3.6
    elif (L/2 + d <= x <= L):
        return 60

def que_const(x, Qp0=400, Qm0=200):
    return Qp0-Qm0

def que(x, Qp0=400, Qm0=200, sigma=0.1, theta=1, L=1):
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
    #print(f"n = {n}")
    h = L/(n+1) #passo dos apoios
    #print(f"h = {h}")

    apoios = [i*h for i in range(n+2)]
    #print("apoios não normalizados: ")
    #print(apoios)

    apoios_barr = [apoios[i]/L for i in range(n+2)]
    #print("apoios normalizados: ")
    #print(apoios_barr)
    #print("\n")
    #print("\n")


    #definição da matriz normal do MMQ 

    #diagonal inferior
    #print("Diagonal inferior")
    a = [0 if i == 0 else (-1)*(1/(L*(h**2)))*gauss(k, apoios_barr[i], apoios_barr[i+1], 2, *x_w(2)) for i in range(n)] 
    #print(a)
    #print("\n")
    a = np.array(a, dtype=np.double)

    #diagonal principal
    #print("Diagonal principal")
    b = [(1/(L*(h**2)))*(gauss(k, apoios_barr[i-1], apoios_barr[i], 2, *x_w(2)) + gauss(k, apoios_barr[i], apoios_barr[i+1], 2, *x_w(2))) for i in range(1, n+1)] 
    #print(b)
    #print("\n")
    b = np.array(b, dtype=np.double)

    #diagonal superior 
    #print("Diagonal superior")
    c = [0 if i == n-1 else (-1)*(1/(L*(h**2)))*gauss(k, apoios_barr[i], apoios_barr[i+1], 2, *x_w(2)) for i in range(n)]
    #print(c)
    #print("\n")
    c = np.array(c, dtype=np.double)

    #array de termos independentes
    #print("Termos independentes")
    d = [(gauss((lambda nos: f(nos)*phi(nos, apoios, h, i)), apoios[i-1], apoios[i], 2, *x_w(2)) + gauss((lambda nos: f(nos)*phi(nos, apoios, h, i)), apoios[i], apoios[i+1], 2, *x_w(2))) for i in range(1, n+1)]
    #print(d)
    #print("\n")
    
    #resolução do sistema linear
    #print("alphas")
    alphas = solve_lin_sys_tridiag(*A_to_LU_tridig(a, b, c), d)
    #print(alphas)
    #print("\n")

    eixox = np.linspace(0, L, 100)
    #print("Vetor solução normalizado:")
    v = [0]*len(eixox)
    for j in range(len(eixox)):
        v[j] = calc_v(alphas, eixox[j], apoios, h, n)
    #print(v)
    #print("\n")
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
    #print("Vetor solução particularizado:")
    #print(u)
    #print("\n")
    return u #vetor solução individualizado

def erro_max(u_teor, u_fem): # devolve o erro maximo entre o valor encontrado pelo fem e o teórico
    erro = u_teor - u_fem
    erro_m = np.max(erro)
    #print(f"Erro máximo = {erro_m}")
    return erro_m

def material_varia(q, L, n): #resolução problemática do exercício 4.4
    return solve_generic(n, L, k_mix, q, 0, 0)

def main():
    n = [7, 15, 31, 63, 127]
    u = [0]*len(n)
    u_teor = gabex()
    for i in range(len(n)):
        u[i] = solve_generic(n[i], 1, kex, qex, 0, 0)
    print("U teórico: ")
    print(u_teor)

    print("U obtido")
    for sol in u:
        print(sol)


if __name__ == "__main__":
    main()