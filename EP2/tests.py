from EP2 import *
from utils import *


if __name__ == "__main__":
    print("EP2 - FÓRMULAS DE INTEGRAÇÃO NUMÉRICA DE GAUSS\n\nARTHUR PEDROSO PORTO BELLI NUSP: 11804608\nLETÍCIA HARUMI FURUSAWA    NUSP: 11965585\n\n")


    print("--------- Exercícios do enunciado ---------")
    print("\n1.1 Volume do cubo de aresta unitária")
    print(f"Para n == 6:\t integral == {gauss_2(f_cubo, 0, 1, c0, d1, 6, *x_w(6))}")
    print(f"Para n == 8:\t integral == {gauss_2(f_cubo, 0, 1, c0, d1, 8, *x_w(8))}")
    print(f"Para n == 10:\t integral == {gauss_2(f_cubo, 0, 1, c0, d1, 10, *x_w(10))}")

    print("\n1.2 Volume do tetraedro de arestas unitárias")
    print(f"Para n == 6:\t integral == {gauss_2(f_tetra, 0, 1, c0, d_tetra, 6, *x_w(6))}")
    print(f"Para n == 8:\t integral == {gauss_2(f_tetra, 0, 1, c0, d_tetra, 8, *x_w(8))}")
    print(f"Para n == 10:\t integral == {gauss_2(f_tetra, 0, 1, c0, d_tetra, 10, *x_w(10))}")

    print("\n2.1 Área do primeiro quadrante [0,1] x [0, 1-x^2]")
    print(f"Para n == 6:\t integral == {gauss_2(one, 0, 1, c0, d_3, 6, *x_w(6))}")
    print(f"Para n == 8:\t integral == {gauss_2(one, 0, 1, c0, d_3, 8, *x_w(8))}")
    print(f"Para n == 10:\t integral == {gauss_2(one, 0, 1, c0, d_3, 10, *x_w(10))}")

    print("\n2.2 Área do primeiro quadrante [0,1] x [0, sqrt(1-y)]")
    print(f"Para n == 6:\t integral == {gauss_2(one, 0, 1, c0, d_32, 6, *x_w(6))}")
    print(f"Para n == 8:\t integral == {gauss_2(one, 0, 1, c0, d_32, 8, *x_w(8))}")
    print(f"Para n == 10:\t integral == {gauss_2(one, 0, 1, c0, d_32, 10, *x_w(10))}")

    print("\n3.1 Volume de e^(y/x) para [0.1, 0.5] x [x^3, x^2]")
    print(f"Para n == 6:\t integral == {gauss_2(f_exp, 0.1, 0.5, c4, d4, 6, *x_w(6))}")
    print(f"Para n == 8:\t integral == {gauss_2(f_exp, 0.1, 0.5, c4, d4, 8, *x_w(8))}")
    print(f"Para n == 10:\t integral == {gauss_2(f_exp, 0.1, 0.5, c4, d4, 10, *x_w(10))}")

    print("\n3.2 Área de e^(y/x) para [0.1, 0.5] x [x^3, x^2]")
    print(f"Para n == 6:\t integral == {gauss_2(f_exp_A, 0.1, 0.5, c4, d4, 6, *x_w(6))}")
    print(f"Para n == 8:\t integral == {gauss_2(f_exp_A, 0.1, 0.5, c4, d4, 8, *x_w(8))}")
    print(f"Para n == 10:\t integral == {gauss_2(f_exp_A, 0.1, 0.5, c4, d4, 10, *x_w(10))}")

    print("\n4.1 Volume da calota esférica de raio 1 e altura 1/4")
    print(f"Para n == 6:\t integral == {gauss_2(f_calota, 0.75, 1, c0, d_calota, 6, *x_w(6))}")
    print(f"Para n == 8:\t integral == {gauss_2(f_calota, 0.75, 1, c0, d_calota, 8, *x_w(8))}")
    print(f"Para n == 10:\t integral == {gauss_2(f_calota, 0.75, 1, c0, d_calota, 10, *x_w(10))}")

    print("\n4.2 Volume de [0, e^(-y^2)] x [-1,1]")
    print(f"Para n == 6:\t integral == {gauss_2(f_error, -1, 1, c0, d_error, 6, *x_w(6))}")
    print(f"Para n == 8:\t integral == {gauss_2(f_error, -1, 1, c0, d_error, 8, *x_w(8))}")
    print(f"Para n == 10:\t integral == {gauss_2(f_error, -1, 1, c0, d_error, 10, *x_w(10))}")