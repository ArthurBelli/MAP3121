# EP1 - MAP3121 - Decomposição LU para matrizes Tridiagonais

O presente programa foi desenvolvido inteira e exclusivamente por Arthur Pedroso Porto Belli e Letícia Harumi Furusawa para o primeiro oferecimento da disciplina MAP3121 - Métodos Numéricos e Aplicações de 2022.

## Requisitos
Seu funcionamento depende de Python 3.9.7 e da biblioteca Numpy 1.12.1.

## Funcionamento
O usuário é apresentado à seguinte mensagem no terminal ao executar o código:

```powershell
>python3 EP1.py
EP1 - MAP3121
    ARTHUR PEDROSO PORTO BELLI NUSP:11804608
    LETÍCIA HARUMI FURUSAWA NUSP:11965585


O CARREGAMENTO DOS DADOS PARA ESSE PROGRAMA DEVE SER FEITO EXCLUSIVAMENTE VIA ARQUIVO EXTERNO '.TXT' QUE DEVE NECESSARIAMENTE ESTAR LOCALIZADO NA MESMA PASTA QUE ESSE SCRIPT.
TAL AQUIVO DEVE SEGUIR A SEGUINTE FORMATAÇÃO:
LINHA 1: VALORES DA DIAGONAL A
LINHA 2: VALORES DA DIAGONAL B
LINHA 3: VALORES DA DIAGONAL C
LINHA 4: VALORES DO VETOR D

Digite o nome do arquivo (com extenção):

```
A seguir, um exemplo de formatação do arquivo a ser carregado com valores adequados:

```powershell
>cat teste.txt
1 1 5 8
2 1 7 3
4 3 9 1
1 1 1 1 
```

Obedecidos os requerimentos de funcionamento, o programa deve imprimir no terminar uma mensagem como à seguir:


```powershell
O sistema foi resolvido pelos seguintes valores:
x_1 = -0.09375
x_2 = 0.36718749999999994
x_3 = 0.2421875
x_4 = -0.28124999999999994
```
*Não foram implementadas medidas de verificação contra *inputs* inadequados. Espera-se que o usuário faça o uso correto do programa.
