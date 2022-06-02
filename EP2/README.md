# EP1 - MAP3121 - Fórmulas de Integração Numérica de Gauss

O presente programa foi desenvolvido inteira e exclusivamente por Arthur Pedroso Porto Belli e Letícia Harumi Furusawa para o primeiro oferecimento da disciplina MAP3121 - Métodos Numéricos e Aplicações de 2022.

## Requisitos
Seu funcionamento depende de Python 3.9.7 e da biblioteca Numpy 1.12.1.

## Funcionamento
O código está separado em três arquivos diferentes:
- ```EP2.py``` é o arquivo principal, que recebe *inputs* do usuário para cálculo de integrais iteradas de funções suaves;
- ```utils.py``` contém funções secundárias;
- ```tests.py``` é o script que contém a resolução dos exercícios pedidos no enunciado. Ele não recebe *inputs*.



O script ```EP2.py``` recebe como *inputs* do usuário: os limites de integração externos ($[a, b]$) e internos ($[c(x), d(x)]$ ou $[c(y), d(y)]$), o número de nós e pesos que se deseja usar ($n$), assim como a função que se deseja integrar $f(x,y)$. Por conta de limitações de tempo e de aplicabilidade, é fornecida uma coleção de classes de funções das quais o usuário pode escolher a que desejar, configurar seus coeficientes e realizar a integral.

$$I = \int_a^b \int_{c(x)}^{d(x)} f(x, y)\space dydx$$
$$I = \int_a^b \int_{c(y)}^{d(y)} f(x, y)\space dxdy$$
