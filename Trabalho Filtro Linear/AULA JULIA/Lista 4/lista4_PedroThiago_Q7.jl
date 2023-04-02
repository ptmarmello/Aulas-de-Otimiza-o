using JuMP,Plots, GLPK, LinearAlgebra

A = [2 1;1 2]
C = [4 3]
b = [4;4]

β = [3;4]
η = [1;2]

MaxCount = 0
## simplex_FaseII e todo seu desenvolvimento pelo algoritmo do Simplex Revisado

function simplex_FaseII(β_init , η_init , A , C , b, MaxCount)

    ## Inicializa o que será utilizado
    B = A[:,β_init]
    N = A[:,η_init]
    cB = C[:,β_init]
    cN = C[:,η_init]
    XB = B\b
    b_barra = B\b

    w = cB*inv(B)
    C_reduzido = -w*N + cN

    V_zero = zeros(size(C_reduzido,2))

    β=β_init
    η=η_init
    ## Loop que verifica o Custo reduzido.
    ## Parte principal do Programa

    while all(C_reduzido .> V_zero)
        ## encontra o menor valor para i (que neste caso foi nomeado k)
        for i in 1:size(N,1)
            a = N[:,i]
            p = w*a-cN[:,i]
            if i == 1
                global save = p
                global k = i
            elseif p < save
                global save = p
                global k = i
            end
        end
        global save = 0
        ## N[:,k] é o que entra na base
        y = inv(B)*N[:,k]

        ## encontra o maior valor para j
        for i in 1:size(y,1)
            if (y[i]>0)
                p = b_barra[i]/y[i]
                if b_barra[i]/y[i] == 0
                    global save = p
                    global l = i
                elseif p > save
                    global save = p
                    global l = i
                end
            end
        end

        ## realiza as devidas mudanças nos índices

        aux = η[k]
        η[k] = size(A,1) + l
        β[l] = aux

        MatrizB = [inv(B) y]
        sizeMB = size(MatrizB,1)

        ultima_col = MatrizB[:,size(MatrizB,2)]
        escolhido = ultima_col[l]

        ## realiza as mudanças na matriz [B^-1 y]
        if escolhido != 0
             escolhido = MatrizB[l,:]/escolhido
        end

        for d in 1:sizeMB
            ## sim aqui dava pra ter melhorado o código
            if d != l
                if (ultima_col[d] != 0)
                    norma = (- ultima_col[d])
                    MatrizB[d,:] = MatrizB[d,:] + norma*MatrizB[l,:]
                    b_barra[d] = b_barra[d] + norma*b_barra[l]
                end
            end
        end

        ## Reestabelece as matrizes que foram utilizadas no passo anterior
        y = MatrizB[:,size(MatrizB,2)]
        new_β= zeros(Int8,size(MatrizB,2)-1)
        for i in 1:size(MatrizB,2)-1
            new_β[i] = i
        end
        B = MatrizB[:,new_β]
        N = A[:,η]
        cB = C[:,β]
        cN = C[:,η]

        ## contador para o número de vezes que o loop acontece
        ## Etapa final:
        w = cB*inv(B)
        C_reduzido = -w*N + cN

        # Devido a troca do auxiliar no passo de troca de índices, teremos essa mudança aqui também.
        β_init = η
        η_init = β
        MaxCount = MaxCount +1
    end # while

    println("Ótimo Encontrado!")
    Xβ= B\b_barra
    Xη=zeros(size(A,2)-size(Xβ,1))
    X=[Xβ;Xη]
    Z = C*X

    return MaxCount, Z, η_init, β_init
end

##

function linear_programming(β, η, A, C, b, MaxCount)


    ## Inicializa as matrízes que serão enviadas ao simplex_FaseII

    Id = diagm(0=>fill(1., size(A,1)))
    A = [A Id]
    Id_C = transpose(fill(0,size(Id,1)))
    C = [C Id_C]

    count, Z, β, η = simplex_FaseII(β, η, A, C, b, MaxCount)

    # Constroi as matrízes utilizando os dados retornados do simplex
    B = A[:,β]
    BT = transpose(B)
    N = A[:,η]
    NT = transpose(N)
    cB = C[:,β]
    cN = C[:,η]
    cTB = transpose(cB)
    cTN = transpose(cN)
    CT = transpose(C)

    X = [1:size(A,1)]
    w = cB*inv(B)
    C_reduzido = w*N - cN

    ## encontra o j máximo para construir dB e constroi dN
    j=findmax(C_reduzido)[2][2]
    dB= -inv(B)*N[:,j]
    dN= -inv(N)*B*dB

    ## Cria as entregas d, X, Z
    d = [dB;dN]

    Xβ=B\b
    Xη=zeros(size(A,2)-size(Xβ,1))
    X=[Xβ;Xη]

    return count, d, X, β, η, Z, C_reduzido
end

##

## Fase 2 para Problema de Produção:

# Seja o problema de produção para:
    #max {CT*X| A*X<=b, X>=0}

count, d, X, β_fim , η_fim, Z, C_reduzido = linear_programming(β, η, A, C, b, MaxCount)




## Letra b
# ##
# θ=zeros(19*19*200+1)
# k = 1
# for j in 2:20
#     global m = j
#     for p in 2:20
#         global n = p
#         for i in 1:200
#             A = rand(m,n)
#             C = rand(1,n)
#             b = round.(rand(m,1))
#
#             η_init = zeros(Int8,n)
#             β_init = zeros(Int8,m)
#
#             for i in 1:n
#                 η_init[i] = Int(i)
#             end
#             for i in 1:m
#                 β_init[i] = Int(n+i)
#             end
#
#             Id = diagm(0=>fill(1., size(A,1)))
#             A = [A Id]
#             Id_C = transpose(fill(0,size(Id,1)))
#             C = [C Id_C]
#
#             Z, β_init, η_init, θ[k] = simplex_FaseII(β_init , η_init , A , C , b)
#             global k+=1
#         end
#     end
# end
#
# ## houve um erro durante o desenvolvimento final para obter os valores do Quantis e por isso não consegui plotar os gráficos
