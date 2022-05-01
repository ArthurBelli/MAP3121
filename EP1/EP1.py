import numpy as np
from math import *

#exercício 1
def LU_to_A(L, U):
    dim = L.shape[0]
    A = np.zeros((dim, dim), dtype=np.double) #matriz A tem mesma dimensão que as de entrada, mas inicializada com zeros
    for i in range(dim): #iterando as linhas
        for j in range(dim): #iterando as colunas
            for k in range(min(i, j)+1): #iterando na regra de formação
                A[i][j] += L[i][k]*U[k][j]
    return A

#exercício 2
def A_to_LU(A):
    dim = A.shape[0]
    L = np.eye(dim, dtype=np.double) #matriz identidade
    U = np.zeros((dim, dim), dtype=np.double)

    for i in range(dim):
        U[i, i:] = A[i, i:] - np.dot(L[i, :i], U[:i, i:])
        L[(i+1):, i] = (A[(i+1):, i] - np.dot(L[(i+1):, :], U[:, i]))/U[i,i]
    return L, U

#exercício 3
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

#retorna as diagonais em forma de vetor de forma mais generalizada, funciona tanto para tridiagonal, como cíclica
def A_to_abc(A):
    dim = A.shape[0]
    a = np.zeros(dim, dtype=np.double)
    b = np.zeros(dim, dtype=np.double)
    c = np.zeros(dim, dtype=np.double)

    for i in range(dim):
        for j in range(dim):
            if i == j:
                b[i] = A[i][j]
            elif i - j == 1 or (i == 0 and j == dim-1):
                a[i] = A[i][j]
            elif i - j == -1 or (i == dim-1 and j == 0):
                c[i] = A[i][j]
    return a, b, c

#retorna a matriz montada a partir de suas diagonais cíclicas
def abc_to_A(a, b, c):
    dim = b.shape[0]
    A = np.zeros((dim, dim), dtype=np.double)
    for i in range(dim):
        for j in range(dim):
            if i == j:
                A[i][j] = b[i]
            elif i - j == 1 or (i == 0 and j == dim-1):
                A[i][j] = a[i]
            elif i - j == -1 or (i == dim-1 and j == 0):
                A[i][j] = c[i]
    return A

#retorna a submatriz diagonal principal de A
def get_T(A):
    dim = A.shape[0]-1 #submatriz diagonal principal
    T = np.zeros((dim, dim), dtype=np.double)
    for i in range(dim):
        for j in range(dim):
            T[i][j] = A[i][j]
    return T

#a pendender do valor de option, retorna o vetor v, ou o vetor w da matrix A
def get_w_v(A, option):
    dim = A.shape[0]-1
    answer = np.zeros(dim, dtype=np.double)
    if option == "v":
        for i in range(dim):
            answer[i] = A[i][dim]
    else:
        for i in range(dim):
            answer[i] = A[dim][i]
    return answer.T #retornamos transposto pois desejamos que seja um vetor coluna

#recebe a matriz A e o vetor d, retorna os valores que resolvem Ax = d
def solve_lin_sys_tridig_cyclic(d, a, b, c):
    A = abc_to_A(a, b, c)
    dim = A.shape[0]
    
    x = np.zeros((dim, 1), dtype=np.double)
    x_barr = np.zeros(dim, dtype=np.double)

    T = get_T(A)
    v = get_w_v(A, "v")
    w = get_w_v(A, "w")

    #T*~y = ~d
    y_barr = solve_lin_sys_tridiag(*A_to_LU_tridig(*A_to_abc(T)), d[0:-1]) #pegamos os valores de d até o penúltimo
    #T*~z = v
    z_barr = solve_lin_sys_tridiag(*A_to_LU_tridig(*A_to_abc(T)), v)

    
    x_n = (d[-1] - c[-1]*y_barr[0]-a[-1]*y_barr[-1])/(b[-1]-c[-1]*z_barr[0]-a[-1]*z_barr[-1])
    x_barr = y_barr - x_n*z_barr

    x = np.append(x_barr, np.array(x_n))
    return x

def read_abc_from_file(path):
    a = []
    b = []
    c = []
    d = []
    with open(path, mode='r', encoding='utf-8') as file:
        line = file.readline().split(" ")
        for value in line: #primeira linha
            a.append(np.double(value))
        line = file.readline().split(" ")
        for value in line: #segunda linha
            b.append(np.double(value))
        line = file.readline().split(" ")
        for value in line: #terceira linha
            c.append(np.double(value))
        line = file.readline().split(" ")
        for value in line: #quarta linha
            d.append(np.double(value))
    return np.array(a, dtype=np.double), np.array(b, dtype=np.double), np.array(c, dtype=np.double), np.array(d, dtype=np.double)

def pretty_print(x):
    print("O sistema foi resolvido pelos seguintes valores: ")
    n = x.shape[0]
    for i in range(n):
        print(f"x_{i+1} = {x[i]}")

def main():
    print("EP1 - MAP3121\n\tARTHUR PEDROSO PORTO BELLI NUSP:11804608\n\tLETÍCIA HARUMI FURUSAWA NUSP:11965585\n\n")
    print("O CARREGAMENTO DOS DADOS PARA ESSE PROGRAMA DEVE SER FEITO EXCLUSIVAMENTE VIA ARQUIVO EXTERNO .TXT QUE DEVE NECESSARIAMENTE ESTAR LOCALIZADO NA MESMA PASTA QUE ESSE SCRIPT.\nTAL AQUIVO DEVE SEGUIR A SEGUINTE FORMATAÇÃO:\nLINHA 1: VALORES DA DIAGONAL A\nLINHA 2: VALORES DA DIAGONAL B\nLINHA 3: VALORES DA DIAGONAL C\nLINHA 4: VALORES DO VETOR D\n")
    path = input("\nDigite o nome do arquivo (com extenção): ")
    a, b, c, d = read_abc_from_file(path)
    x = solve_lin_sys_tridig_cyclic(d, a, b, c)
    pretty_print(x)

if __name__ == "__main__":
    main()