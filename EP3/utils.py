'''
Como dependemos de funções implementadas nos outros EPs, 
as funções necessárias foram copiadas e, a depender do caso, adaptadas para esse arquivo a 
fim de facilitar a importação das funcionalidades sem a necessidade 
de se ter os arquivos originais no mesmo diretório
'''

# Funções do EP1
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

#retorna a matriz montada a partir de suas diagonais
def abc_to_A(a, b, c):
    dim = b.shape[0]
    A = np.zeros((dim, dim), dtype=np.double)
    for i in range(dim):
        for j in range(dim):
            if i == j:
                A[i][j] = b[i]
            elif i - j == 1:
                A[i][j] = a[i]
            elif i - j == -1:
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

#retorna valor v da matriz passada
def get_v(A, option):
    dim = A.shape[0]-1
    answer = np.zeros(dim, dtype=np.double)
    for i in range(dim):
        answer[i] = A[i][dim]
    return answer.T #retornamos transposto pois desejamos que seja um vetor coluna


#recebe a matriz A e o vetor d, retorna os valores que resolvem Ax = d
def solve_lin_sys_tridig_cyclic(d, a, b, c):
    A = abc_to_A(a, b, c)
    dim = A.shape[0]
    
    x = np.zeros((dim, 1), dtype=np.double)
    x_barr = np.zeros(dim, dtype=np.double)

    T = get_T(A)
    v = get_v(A)

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


# Funções do EP2

## Funções do arquivo utils.py

def c0(x):
    return 0

def c(f, x):
    return f(x)

def c4(x):
    return x**3

def d4(x):
    return x**2

def d_3(x):
    return 1-x**2

def d_32(x):
    return np.sqrt(1-x)
    
def d1(x):
    return 1

def d_calota(x):
    return np.sqrt(1-x**2)

def d_tetra(x):
    return 1-x

def f_cubo(x, y):
    return 1

def f_tetra(x, y):
    return 1-x-y

def f_exp(x, y):
    return np.exp(y/x)

def one(x, y):
    return 1

def f_exp_A(x, y):
    return np.sqrt(1+ ((-1)*(y/x**2)*f_exp(x, y))**2 + (f_exp(x, y)/x)**2)

def f_calota(x, y):
    return 2*np.pi*y

def d_error(x):
    return np.exp((-1)*(x**2))

def f_error(x, y):
    return 2*np.pi*d_error(x)

def exp(x, y=0):
    return np.exp(x)

def exp2(x, y=0):
    return np.exp((x**2))

def expxy(x, y=0):
    return np.exp(y/x)

def log(x, y=0):
    return np.log10(x)

def ln(x, y=0):
    return np.log(x)

def poly(x, y=0): # função polinomial arbitrária para fins de demonstração
    return x**4 + 3*x**2 - x + 12

def cos(x, y=0): #recebe x em radianos
    return np.cos(x)

def sin(x, y=0):
    return np.sin(x)

def one(x, y=0):
    return 1

def zero(x, y=0):
    return 0

funcas = { "e^x": exp, "e^(x^2)": exp2, "e^(y/x)": expxy, "log(x)": log, "ln(x)": ln, "polynomial(x)": poly, "cos(x)": cos, "sin(x)": sin, "one": one, "zero": zero,
           "1-x-y": f_tetra, "1-x": d_tetra, "x^2": d4, "x^3": c4, "1-x^2": d_3, "sqrt(1-x)": d_32, "sqrt(1-x^2)": d_calota, }


## Funções do arquivo EP2.py

def x_w(n):
    if n == 6:
        x = [-0.9324695142031520278123016,-0.6612093864662645136613996,-0.2386191860831969086305017,
             0.2386191860831969086305017,0.6612093864662645136613996,0.9324695142031520278123016]
        w = [0.1713244923791703450402961,0.3607615730481386075698335,0.4679139345726910473898703,
             0.4679139345726910473898703,0.3607615730481386075698335,0.1713244923791703450402961]
    
    if n == 8 :
        x = [-0.9602898564975362316835609,-0.7966664774136267395915539,-0.5255324099163289858177390,
             -0.1834346424956498049394761,0.1834346424956498049394761,0.5255324099163289858177390,
             0.7966664774136267395915539,0.9602898564975362316835609]
        w = [0.1012285362903762591525314,0.2223810344533744705443560,0.3137066458778872873379622,
             0.3626837833783619829651504,0.3626837833783619829651504,0.3137066458778872873379622,
             0.2223810344533744705443560,0.1012285362903762591525314]

    if n == 10:
        x = [-0.9739065285171717200779640,-0.8650633666889845107320967,-0.6794095682990244062343274,
             -0.4333953941292471907992659,-0.1488743389816312108848260,0.1488743389816312108848260,
             0.4333953941292471907992659,0.6794095682990244062343274,0.8650633666889845107320967,0.9739065285171717200779640]
        w = [0.0666713443086881375935688,0.1494513491505805931457763,0.2190863625159820439955349, 
             0.2692667193099963550912269,0.2955242247147528701738930,0.2955242247147528701738930,
             0.2692667193099963550912269,0.2190863625159820439955349,0.1494513491505805931457763,0.0666713443086881375935688]
    
    return x, w

# formula da gauss com mudança de variável


def find_lin_scale_transp(a, b):
    '''
    função retorna os fatores de correção de intervalos [a,b] para [-1,1]
    '''
    return np.double((b-a)/2), np.double((a+b)/2)

def gauss(f, a, b, n, x, w):
    '''
    integração simples no intervalo [a,b] da função f(x) com n nós (x) e pesos (w)
    '''
    linear_scaling, linear_transposition = find_lin_scale_transp(a, b)
    sum = 0
    v = [linear_scaling*w[i] for i in range(n)] # pesos adaptados
    y = [(linear_scaling*x[i]+linear_transposition) for i in range(n)] # nós adaptados
    for i in range (n): # soma iterada da integral
        sum += v[i]*f(y[i])

    return sum


#foi necessário alterar a função, modo que agora ela recebe a função envelope f, a função a ser passada para f (g) e o modo de operação, 
#que determina se deve ser utilizado x_{i-1} ou x_{i+1}

def gauss_2(f, g, a, b, c, d, n, nos, pesos, mode):
    '''
    integração dupla no intervalo [a,b] x [c(x_i), d(x_i)] da função f(x, y) com n nós (x) e pesos (w).
    No caso do EP3, nossas funções são unidimensionais, ou seja, dependem apenas de x
    '''
    linear_scaling_ext, linear_transposition_ext = find_lin_scale_transp(a, b) # fatores de correção para a integral externa
    pesos_ext = [linear_scaling_ext*pesos[i] for i in range(n)] # pesos da integral externa
    nos_ext = [(linear_scaling_ext*nos[i]+linear_transposition_ext) for i in range(n)] # nós da integração externa

    sum_ext = 0
    for i in range(n): #para cada nó externo, definiremos o intervalo de integração, seus nós e seus pesos
        a_int = c(nos_ext[i]) # extremo inferior de integração interno a cada iteração
        b_int = d(nos_ext[i]) # extremo superior de integração interno a cada iteração

        linear_scaling_int, linear_transposition_int = find_lin_scale_transp(a_int, b_int) # fatores de correção para a integral interna a cada iteração
        pesos_int = [linear_scaling_int*pesos[i] for i in range(n)] # pesos da integral interna
        nos_int = [linear_scaling_int*nos[i]+linear_transposition_int for i in range(n)] # nós da integral interna
        
        sum_int = 0
        for j in range(n):
            if mode == "inf":
                sum_int += pesos_int[j]*f(g, nos_ext[i], nos_int[j], c) #principal alteração de funcionalidade
            else:
                sum_int += pesos_int[j]*f(g, nos_ext[i], nos_int[j], d)
        sum_ext += pesos_ext[i]*sum_int
    
    return sum_ext
