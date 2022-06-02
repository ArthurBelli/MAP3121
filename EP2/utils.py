import numpy as np

'''
Funções para cálculo dos exercícios propostos no enunciado
'''

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
    return np.exp((-1)*x**2)

def f_error(x, y):
    return 2*np.pi*abs(x)

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

#função polinomial para expoentes >= 1
def poly(x, y=0): #coefs recebe array de coeficientes do menor para o maior [a_0, a_1, a_2, ..., a_n]
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

def print_funcas():
    print(*funcas, sep='\n')

