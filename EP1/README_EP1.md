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


## Execução
O usuário será aprentado com a seguinte sequência de mensagens ao executar o ```EP2.py```:
```powershell
>python .\EP2.py
e^x
e^(x^2)
e^(y/x)
log(x)
ln(x)
polynomial(x)
cos(x)
sin(x)
one
zero
1-x-y
1-x
x^2
x^3
1-x^2
sqrt(1-x)
sqrt(1-x^2)
Escolha a função da lista que deseja integrar:
```

Para fazer a seleção da função, basta escrever a expressão desejada.

Passado disso, serão pedidos os extremos de integração internos e em sequência os extremos externos. Após isso, será apresentado o resultado final para os 3 valores possíveis de $n$: 6, 8 e 10 nós.

A título de exemplo, será feita a integral: $$\int_5^{10} \int_0^{e^x} e^x\space dydx$$

```powershell
> python .\EP2.py
e^x
e^(x^2)
e^(y/x)
log(x)
ln(x)
polynomial(x)
cos(x)
sin(x)
one
zero
1-x-y
1-x
x^2
x^3
1-x^2
sqrt(1-x)
sqrt(1-x^2)
Escolha a função da lista que deseja integrar: e^x
Da mesma lista de funções, escolha agora o extremo de integração superior: e^x
Da mesma lista de funções, escolha agora o extremo de integração inferior: zero
O extremo de integração superior exterior: 10
O extremo de integração inferior exterior: 5
n == 6
int_(5.0)^(10.0) int_(zero)^(e^x) e^x dydx = 242566823.26707447
n == 8
int_(5.0)^(10.0) int_(zero)^(e^x) e^x dydx = 242571580.56276858
n == 10
int_(5.0)^(10.0) int_(zero)^(e^x) e^x dydx = 242571584.47075567
```

O resultado obtido é condizente com o resultado exato. Utilizando a calculadora online *Wolfram Alpha* para validar a resposta, encontramos que $$\int_5^{10} \int_0^{e^x} e^x\space dydx = 2.4257158447199773\cdot 10^8$$
