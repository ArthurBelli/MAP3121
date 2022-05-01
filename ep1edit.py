import numpy as np
import math


def LU_to_A(L, U, dim):
    
    A = np.zeros((dim, dim), dtype=np.double) #matriz A tem mesma dimensão que as de entrada, mas inicializada com zeros
    for i in range(dim): #iterando as linhas
        for j in range(dim): #iterando as colunas
            for k in range(min(i, j)+1): #iterando na regra de formação
                A[i][j] += L[i][k]*U[k][j]
    return A

def A_to_LU(A,dim):
    
    L = np.eye(dim, dtype=np.double) #matriz identidade
    U = np.zeros((dim, dim), dtype=np.double)

    for i in range(dim):
        U[i, i:] = A[i, i:] - np.dot(L[i, :i], U[:i, i:])
        L[(i+1):, i] = (A[(i+1):, i] - np.dot(L[(i+1):, :], U[:, i]))/U[i,i]
    
    print(U)
    print()
    print(L)

    L.flatten()
    U.flatten()

    return L, U

def get_abc(A,dim):

    a = np.zeros(dim,dtype=np.double)
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

def solve_lin_sys_tridig(L,U,A,d,dim):

    y = np.zeros(dim, dtype=np.double)
    x = np.zeros(dim, dtype=np.double)

    y[0]=d[0]

    a,b,c = get_abc(A,dim)

    for i in range(1,dim):
       y[i] = d[i]-L[i]*y[i-1]

    x[dim-1] = y[dim-1]/U[dim-1]

    for i in range(dim-2,-1,-1):
       x[i]=(y[i]-c[i]*x[i+1])/U[i]
    
    return(x)

def get_v(A):
    dim = A.shape[0]-1
    v = np.zeros(dim, dtype=np.double)

    for i in range(dim):
        v[i] = A[dim][i]

def solve_lin_sys_cyclic_tridig(L,U,A,d,dim):

    x = np.zeros(dim,dtype=np.double)
    y = np.zeros(dim,dtype=np.double)

    a,b,c = get_abc(A,dim)


    v = get_v(A)

    y_barr = solve_lin_sys_tridig(L,U,A,d,dim-1)
    z_barr = solve_lin_sys_tridig(L,U,A,v,dim-1)

    x[dim-1]=(d[dim-1]-c[dim-1]*y_barr[0]-a[dim-1]*y_barr[dim-2])/(b[dim-1]-c[dim-1]*z_barr[0]-a[dim-1]*z_barr[dim-2])

    for i in range(dim-2,-1,-1):
       x[i]=y_barr[i]-x[-1]*z_barr[i]

    return x

def main() -> None:

    A = np.array([[1, 1, 0, 0, 6], [1, 3, 2, 0, 0], [0, 1, 5, 9, 0], [0, 0, 15, 4, 10], [7, 0, 0, 2, 1]])
    d=([1,1,1,1,1])
    dim=5

    L,U = A_to_LU(A,dim)

    solve_lin_sys_tridig(L,U,A,d,dim)
    
    




    
  
if __name__ == '__main__':
    main()


