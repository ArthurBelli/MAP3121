from ep3 import *
import matplotlib.pyplot as plt

print("MAP3121: EP3 - Modelo de Resfriamento de Chips")
print("ARTHUR PEDROSO PORTO BELLI NUSP 11804608\nLETÍCIA HARUMI FURUSAWA NUSP 11965585")

print("--------- Exercícios do enunciado ---------")
print("1. f(x) = 12x(1-x)-2, k=1, n = 7, L = 1")

n = 7
L = 1
eixox = np.linspace(0, L, 100)
u_teorico = [gabarito(eixox[k]) for k in range(len(eixox))]
u_obtido_7 = solve_generic(7, L, k, excitacao, 0, 0)
plt.plot(eixox, u_teorico, 'r')
plt.plot(eixox, u_obtido_7, 'b')
plt.title(f"Exercício 4.2 para k = 1, n = {n} e L = {L}")
plt.xlabel("L")
plt.ylabel("Temperatura")
plt.legend(["u_teórico", "u_obtido_7"], loc="upper left")
plt.savefig(f"4.2_n{n}L{L}.png")

plt.clf()

print("2. f(x) = 12x(1-x)-2, k = 1, n = 15, L = 2")

n = 15
L = 2
eixox = np.linspace(0, L, 100)
u_obtido_15 = solve_generic(n, L, k, excitacao, 0, 0)
plt.plot(eixox, u_obtido_15, 'b')
plt.title(f"Exercício 4.2 para k = 1, n = {n} e L = {L}")
plt.xlabel("L")
plt.ylabel("Temperatura")
plt.legend(["u_obtido_15"], loc="upper left")
plt.savefig(f"4.2_n{n}L{L}.png")

plt.clf()

print("3. f(x) = 12x(1-x)-2, k = 1, n = 31, L = 1")

n = 31
L = 1
eixox = np.linspace(0, L, 100)
u_teorico = [gabarito(eixox[k]) for k in range(len(eixox))]
u_obtido_31 = solve_generic(n, L, k, excitacao, 0, 0)
plt.plot(eixox, u_teorico, 'r')
plt.plot(eixox, u_obtido_31, 'b')
plt.title(f"Exercício 4.2 para k = 3.6, n = {n} e L = {L}")
plt.xlabel("L")
plt.ylabel("Temperatura")
plt.legend(["u_teorico", "u_obtido_31"], loc="upper left")
plt.savefig(f"4.2_n{n}L{L}.png")

plt.clf()

print("4. f(x) = 12x(1-x)-2, k = 3.6, n = 63, L = 3")

n = 63
L = 3
eixox = np.linspace(0, L, 100)
u_obtido_63 = solve_generic(n, L, k36, excitacao, 0, 0)
plt.plot(eixox, u_obtido_63, 'b')
plt.title(f"Exercício 4.2 para k = 3.6, n = {n} e L = {L}")
plt.xlabel("L")
plt.ylabel("Temperatura")
plt.legend(["u_obtido_63"], loc="upper left")
plt.savefig(f"4.2_n{n}L{L}.png")

plt.clf()

print("5. Exercício com forçantes de calor")
print("Parâmetros: Qp0 = 600, Qm0 = 300, sigma = 0.1, theta = 1, n = 63, L = 1, k = 1")
print("CONDIÇÕES DE CONTORNO NÃO HOMOGÊNEAS: a = 0 e b = 10")

n = 63
L = 1
eixox = np.linspace(0, L, 100)
u_obtido = solve_generic(63, 1, k, que, 1, 10)
plt.plot(eixox, u_obtido, 'b')
plt.title(f"4.3 para Qp0 = 600, Qm0 = 300, sigma = 0.1, theta = 1, n = 63, \nL = 1, k = 1, a = 0, b = 10")
plt.xlabel("L")
plt.ylabel("Temperatura")
plt.legend(["u_obtido"], loc="upper left")
plt.savefig(f"4.3_n{n}L{L}.png")
