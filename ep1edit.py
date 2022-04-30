from typing import List
import numpy as np


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
    
    print(U)
    print()
    print(L)
    return [L, U]

#exercício 3
def A_to_LU_tridig(a, b, c):
    dim = b.shape[0] #número de elementos na diagonal principal é a dimmensão da matriz
    l_ip1 = np.zeros(dim, dtype=np.double) #l_i+1
    u_ii = np.zeros(dim, dtype=np.double) #u_ii
    u_ip1 = c.copy() #u_i+1

    u_ii[0] = b[0]
    for i in range(1, dim):
        l_ip1[i] = a[i]/u_ii[i-1] #multiplicadores
        u_ii[i] = b[i]-l_ip1[i]*c[i-1]

    return u_ii, u_ip1, l_ip1


def solve_lin_sys(u, u_ip1, l, d):
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

def get_abc(A): #retorna as diagonais em forma de vetor
    dim = A.shape[0]
    a = np.zeros(dim, dtype=np.double)
    b = np.zeros(dim, dtype=np.double)
    c = np.zeros(dim, dtype=np.double)

    for i in range(dim):
        for j in range(dim):
            if i == j:
                b[i] = A[i][j]
            elif i - j == 1:
                a[i] = A[i][j]
            elif i - j == -1:
                c[i] = A[i][j]
    return a, b, c

def get_T(A):
    dim = A.shape[0]-1 #submatriz diagonal principal
    T = np.zeros((dim, dim), dtype=np.double)
    for i in range(dim):
        for j in range(dim):
            T[i][j] = A[i][j]
    return T

def get_v(A):
    dim = A.shape[0]-1
    v = np.zeros(dim, dtype=np.double)

    for i in range(dim):
        v[i] = A[dim][i]
    
    return v.T #retornamos a transposta pois desejamos que seja um vetor coluna

def get_w(A):
    dim = A.shape[0]-1
    w = np.zeros(dim, dtype=np.double)

    for i in range(dim):
        w[i] = A[i][dim]

    return w.T #mesmo raciocínio do vetor v

#ta quebrado ainda, mas vai funfar
def solve_lin_sys_tridig(LU:List[np.array],A, d): #"main"
    dim = A.shape[0]
    x = np.zeros(dim, dtype=np.double)
    y = np.zeros(dim, dtype=np.double)

    y[0]= d[0]

    for i in range (dim):
        y[i]= d[i] - LU[0][i]*y[i-1]

    x[-1]=y[-1]/LU[1][-1]
    
    for i in range (dim-2,-1,-1):
        x[i] = (y[i]-A[2][i]*x[i+1]/LU[1][i])

    return x


def solve_lin_sys_cyclic_tridig(A,d):
    dim = A.shape[0]
    x = np.zeros(dim, dtype=np.double)
    y = np.zeros(dim-1, dtype=np.double)

    T = get_T(A)
    v = get_v(A)
    w = get_w (A)

    LU = A_to_LU(T)

    z_barr = solve_lin_sys_tridig(LU,A,v)

    y_barr = solve_lin_sys_tridig(LU,A,d)

    cn = A[dim][0]
    an = A[0][dim]
    bn = A[dim][dim]

    x[dim] = (d[dim] - (cn*y_barr[1]) - (an*y_barr[dim-1]))/(bn - (cn*z_barr[1])- an*z_barr[dim-1])

    for i in range(dim-1):
        x[i]=y_barr[i]-x[dim]*z_barr[i]
        
    return x 

def main() -> None:

    A = np.array([[3,-0.1,-0.2],[0.1,7,-0.3],[0.3,-0.2,10]])
    A_to_LU(A)
    
if __name__ == '__main__':
    main()
