using CSV, Statistics, RCall

vec = [1,2,3]

vec = [i for i=1:100]

vec_2 = [i for i=1:10:1000]

vec_multiplo_7 = [i for i=1:10000 if i%7 == 0]


vec*vec_2
vec*vec_2'

vec .* vec_2

A = [
    1 2 3;
    3 2 1;
]

B = [
    0 0 0;
    3 2 1;
]

A*vec

## Generates a rand matrix
A = rand(100,100)

#Access first column of A
vec = A[:,1]
A[:,1] = zeros(100)
A[:,2] = ones(100)

μ = [2,5]

A[:,[2,5]]

B = A[:,1]

B = ones(100)

mean(A)
std(A)


lin,col = size(A)

for i=1:lin
    A[i,3] = i
end

## Ou

A[:,3] = [i for i=1:lin]


function modifica_matriz(matriz, μ)
    lin,col = size(matriz)
    n_matriz = matriz
    for indice in μ
        n_matriz[:,indice] = zeros(lin)        
    end
    return n_matriz
end

A = rand(10,10)
B = modifica_matriz(A,[2,3])
A

### CAUTION !! immutable types
A = rand(10,10)
isimmutable(A)

### CAUTION !! immutable types

function modifica_matriz_mantendo_original(matriz)
    lin,col = size(matriz)
    n_matriz = matriz[:,:]
    for indice in μ
        n_matriz[:,indice] = zeros(lin)        
    end
    return n_matriz
end

A = rand(10,10)
B = modifica_matriz_mantendo_original(A)
A


## Using strings
str = "π is about 3.1415"

str[1:10]

str_splited = split(str," ")


## Dictionary
preço = Dict(
    "Arroz" => 4,
    "Arroz Integral" => 5,
    "Kinder" => 8,
)

preço["Arroz"]
preço["Kinder"]

vec_strings = ["$i" for i=1:10]

keys(preço)
values(preço)

## Other stuff
♂
◂
Inf

## Why is julia fast?
soma(x) = x + x

@code_warntype soma(1)

@code_warntype soma(1.1)