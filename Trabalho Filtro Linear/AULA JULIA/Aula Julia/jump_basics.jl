## Toda primeira vez usando um pacote, devemos
## devemos adicionar ele usando o Pkg

## Passo 1:
using Pkg

## Passo 2:
##Adicionar o pacote dejesado
##Pkg.add("PACOTE DESEJADO")

## Por exemplo, pacotes interessantes de já terem acesso
## CSV é um pacote para leitura e escrita de arquivos CSV
Pkg.add("CSV")

## DataFrames é essencial para o uso de Dados com muita
## informação
Pkg.add("DataFrames")

## Statistics é útil para análises estatísticas.. 
## cuidade de a função média não é nativa no julia.. usar 
## o pacote Statistics
Pkg.add("Statistics")

## JuMP é o pacote de otimização. É uma interface
## para escrever modelos de otimização de forma intuitiva.
Pkg.add("JuMP")

## GLPK é um solver de problemas lineares. O modelo 
## criado através do JuMP é ordenado e enviado ao GLPK
## que retorna o resultados do problema.
Pkg.add("GLPK")

## Benchmark é um pacote para avaliar performance de 
## códigos. Retorna tempo médio de realização. Alocaçoes
## memmória utilizada..
Pkg.add("Benchmark")



## Para definir um modelo de otimização, basta declarar:
modelo = Model()

## Porém, devemos passar também o solver que ira resolver
## o modelo:
modelo = Model(with_optimizer = GLPK.optimizer)

## Podemos adicionar váriaveis ao problema, passando sempre
## o modelo que estamos usando
@variable(modelo, x1)
@variable(modelo, x2)

## Cuidade para não atribuir o mesmo nome as variaveis

## Podemos declarar a variavel já com algumas restrições.
## Por exemplo, uma variavel sempre positiva ou negativa
## Ou ainda uma variavel inteira, ou binária. 
## Lembrem que problemas lineares(convexos) apenas admitem variaveis
## contínuas.
@variable(modelo, x1>=0)
@variable(modelo, x1<=0)
@variable(modelo, x1, Bin)
@variable(modelo, x1, Int)

## Para adiocionar restrições use:
@contraint(modelo, 2*x1 + x2 <= 4)
@contraint(modelo, x1 + 4*x2 <= 4)


## Caso precisemos criar restrições/variaveis iterativamente podemos:
## Definir variavel dando um intervalo
@variable(modelo, x[1:10])

## As variaveis do jump, ao serem declaradas ficam associadas ao indice
## que ela foi definida. Isto é podemos criar uma variavel da forma:
@variable(modelo, x[[1,5,10]])

## Apenas 3 variaveis foram criadas, x[1], x[5] e x[10]. E se quisermos 
## acessa-las basta usar esse indice
x[5]

## Isso vale porque as variaveis JuMP são criadas da mesma forma de um dicionario,
## cada indice é como uma chave apontando para a variavel. Por exemplo:
modelo = Model()
@variable(modelo, x[["hello","i","exist"]])



## Voltando a declarar iterativamente
## Restrições dentro de um for/while
for i=1:10
    @contraint(modelo, x[i] <= 4)
end

## Ou ainda
@contraint(modelo, x[1:10]<=4)

## Ou ainda..
@contraint(modelo, x[["hi","i","exist"]]<=4)


## Lembrando que colocar . antes do operador realiza para 
## cada elemento..
@contraint(modelo, x.<=4)


## Função Objetivo é criada, suponha x como um vetor de variaveis
## com 10 elementos.
## Queremos Minimizar ou Maximizar a função objetivo
@objective(modelo, Min, sum(x for i=1:10))
@objective(modelo, Max, sum(x for i=1:10))



## Ou ainda 
@objective(modelo, Min, x[1]+x[2]+x[3]+x[4]))



## Para otimizar o modelo JuMP usamos
optimize!(modelo)

## E acessar o status do resultados
termination_status(modelo)

## Valor da Função objetiva
objective_value(modelo)

## Valor das variaveis.. dependendo da forma que você as criou

value(x1)
value(x[1])
value.(x)


## Para mais sobre como usar o JuMP, acesse 
# https://github.com/JuliaOpt/JuMP.jl -  Documentation

# Material preparado pelo LAMPS PUC Rio
# https://github.com/LAMPSPUC/Teaching.jl/tree/master/Optimization


#### PROBLEMA DA PRODUÇÃO - EXEMPLO

produção = Model(with_optimizer = GLPK.optimizer)

@variable(produção, x1)
@variable(produção, x2)

@constraint(produção, 2*x1 + x2 <= 4)
@constraint(produção, x1 + 2*x2 <= 4)

@objective(produção, Max, 4*x1 + 3*x2)

optimize!(produção)

objective_value(produção)

value(x1)
value(x2)


## Forma Padrão
c = [4 3]
A = [2 1;1 2]
b = [4 4]

produção = Model(with_optimizer = GLPK.optimizer)

@variable(produção, x[1:2])

@constraint(produção, A*x <= b)

@objective(produção, c*x)

optimize!(produção)

objective_value(produção)

value.(x)