import numpy as np
import utils

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



def gauss_2(f, a, b, c, d, n, x, w):
    '''
    integração dupla no intervalo [a,b] x [c(x_i), d(x_i)] da função f(x, y) com n nós (x) e pesos (w)
    '''
    linear_scaling_ext, linear_transposition_ext = find_lin_scale_transp(a, b) # fatores de correção para a integral externa
    pesos_ext = [linear_scaling_ext*w[i] for i in range(n)] # pesos da integral externa
    nos_ext = [(linear_scaling_ext*x[i]+linear_transposition_ext) for i in range(n)] # nós da integração externa

    sum_ext = 0
    for i in range(n): #para cada nó externo, definiremos o intervalo de integração, seus nós e seus pesos
        a_int = c(nos_ext[i]) # extremo inferior de integração interno a cada iteração
        b_int = d(nos_ext[i]) # extremo superior de integração interno a cada iteração

        linear_scaling_int, linear_transposition_int = find_lin_scale_transp(a_int, b_int) # fatores de correção para a integral interna a cada iteração
        pesos_int = [linear_scaling_int*w[i] for i in range(n)] # pesos da integral interna
        nos_int = [linear_scaling_int*x[i]+linear_transposition_int for i in range(n)] # nós da integral interna
        
        sum_int = 0
        for j in range(n):
            sum_int += pesos_int[j]*f(nos_ext[i], nos_int[j])
        sum_ext += pesos_ext[i]*sum_int
    
    return sum_ext

def main():
    
    utils.print_funcas()
    f = input("Escolha a função da lista que deseja integrar: ")
    d = input("Da mesma lista de funções, escolha agora o extremo de integração superior: ")
    c = input("Da mesma lista de funções, escolha agora o extremo de integração inferior: ")
    b = np.double(input("O extremo de integração superior exterior: "))
    a = np.double(input("O extremo de integração inferior exterior: "))

    for n in [6, 8, 10]:
        I = gauss_2(utils.funcas[f], a, b, utils.funcas[c], utils.funcas[d], n, *x_w(n))
        print(f"n == {n}")
        print(f"int_({a})^({b}) int_({c})^({d}) {f} dydx = {I}")

if __name__ == "__main__":
    main()