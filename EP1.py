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
def get_abc(A):
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
def solve_lin_sys_tridig_cyclic(A, d):
    dim = A.shape[0]
    x = np.zeros((dim, 1), dtype=np.double)
    x_barr = np.zeros(dim, dtype=np.double)
    a, b, c = get_abc(A)

    T = get_T(A)
    v = get_w_v(A, "v")
    w = get_w_v(A, "w")

    #T*~y = ~d
    y_barr = solve_lin_sys_tridiag(*A_to_LU_tridig(*get_abc(T)), d[0:-1]) #pegamos os valores de d até o penúltimo
    #T*~z = v
    z_barr = solve_lin_sys_tridiag(*A_to_LU_tridig(*get_abc(T)), v)

    
    x_n = (d[-1] - c[-1]*y_barr[0]-a[-1]*y_barr[-1])/(b[-1]-c[-1]*z_barr[0]-a[-1]*z_barr[-1])
    x_barr = y_barr - x_n*z_barr

    x = np.append(x_barr, np.array(x_n))
    return x

# testes
A = np.array([[2, 4, 0, 7], [2, 7, 4, 0], [0, 7, 9, 5], [1, 0, 3, 8]])
d = np.array([2, -1, 0, 5])


a = [(2*i-1)/(4*i) for i in range(1, 20)]
a.append((2*20-1)/(2*20))
print(a)


print(solve_lin_sys_tridig_cyclic(A, d))