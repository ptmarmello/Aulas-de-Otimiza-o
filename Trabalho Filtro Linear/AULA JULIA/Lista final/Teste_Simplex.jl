using JuMP,Plots, GLPK, LinearAlgebra

A = [2 1;1 2]
C = [4 3]
b = [4;4]

β = [3;4]
η = [1;2]

MaxCount = 0
function simplex_FaseII(β_init , η_init , A , C , b, MaxCount)

    Id = diagm(0=>fill(1., size(A,1)))
    A = [A Id]
    Id_C = transpose(fill(0,size(Id,1)))
    C = [C Id_C]

    ## Inicializa o que será utilizado
    B = A[:,β_init]
    N = A[:,η_init]
    cB = C[:,β_init]
    cN = C[:,η_init]
    XB = B\b
    b_barra = B\b

    w = cB*inv(B)
    C_reduzido = -w*N + cN

    global pass = 0
    for count in 1:size(C_reduzido,2)
        if C_reduzido[count] .<= 0
            global pass +=1
        end
    end
    if pass == size(C_reduzido,2)
        global pass = false
    else
        global pass = true
    end

    β=β_init
    η=η_init
    ## Loop que verifica o Custo reduzido.
    ## Parte principal do Programa
    while pass
        save = 0
        for counter in 1:size(C_reduzido,2)
            temp = C_reduzido[counter]
            if counter == 1
                global j = counter
                global save = temp
            elseif temp > save
                global j = counter
                global save = temp
            end
        end

        a = N[:,j]
        dB_linha = -B\a

        if all(dB_linha .>= zeros(size(dB_linha,2)) )
            return ("ilimitado")
        else

            global save = 0
            for i in 1:size(XB,1)
                temp = XB[i]/abs(dB_linha[i])
                if i == 1
                    global k = 1
                    global save = temp
                elseif temp < save
                    global k = i
                end # if
            end #for
        end #if

        aux = η[j]
        η[j] = β[k]
        β[k] = aux
        B = A[:,β]
        XB = B\b
        N = A[:,η]
        cB = C[:,β]
        cN = C[:,η]
        C_reduzido = -w*N + cN
        MaxCount = MaxCount +1
        global pass = 0
        for count in 1:size(C_reduzido,2)
            if C_reduzido[count] .<= 0
                global pass +=1
            end
        end
        if pass == size(C_reduzido,2)
            global pass = false
        else
            global pass = true
        end

    end # while

    println("Ótimo Encontrado!")
    Xβ= B\b_barra
    Xη=zeros(size(A,2)-size(Xβ,1))
    X=[Xβ;Xη]
    Z = C*X

    return MaxCount, Z, η, β
end

##

function linear_programming(β, η, A, C, b, MaxCount)

    return simplex_FaseII(β, η, A, C, b, MaxCount)

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


## Função principal que organiza a ida para as fases I e II
function simplex_Programming(A , C, b, MaxCount)

    ## Análise::
    V_zero = zeros(size(b))
    n = size(A,2)
    m = size(A,1)


    ## Verifica se b >=0, dessa forma, sendo direcionado para fase II
    if all(b .>= V_zero)
        η_init = zeros(Int8,n)
        β_init = zeros(Int8,m)

        for i in 1:n
            η_init[i] = Int(i)
        end
        for i in 1:m
            β_init[i] = Int(n+i)
        end
        ## União das Matrízes e Preparação para as Fases::
        return simplex_FaseII(β_init, η_init, A, C, b, MaxCount)
    end

    ## Caso vá para a FaseI
    return simplex_FaseI(A,C,b,MaxCount)
end


## Início Fase I
function simplex_FaseI(A,C,b, MaxCount)

    ## n é o número de colunas, m é o número de linhas de A_barra
    n=size(A,2)
    m=size(A,1)

    # Inicializa a matriz e (que multiplica w) e a nova matriz A
    ϵ = fill(-1.,m)
    Ae = [A ϵ]

    #Inicializa o eta e beta que marcarão as posições
    η_init = zeros(n+1)
    β_init = zeros(m)

    for i in 1:n
        η_init[i] = i
    end

            ## Essa parte não precisa, já que vai substituir de qualquer forma
            η_init[n+1] = n+m+1

    for i in 1:m
        β_init[i] = n+i
    end

    # θ=zeros(length(b))
    for i in 2:length(b)
        if b[i]< b[i-1]
            global θ=i
        else
            global θ=i-1
        end
    end

    aux = β_init[θ]
    β_init[θ] = η_init[n+1]
    η_init[n+1] = aux


    w = zeros(n+m+1)
    w[n+m+1]= -1.

    ### Call fase II para w

    Z_faseI, β_faseI, η_faseI = simplex_FaseII( β_init, η_init, Ae, -w, b, MaxCount)
    tol = 10^-5
    # Em algum momento isso aqui deve ser feito:
    if (Z_faseI == 0 || Z_faseI > -tol)
        return simplex_FaseII( β_faseI, η_faseI, A, C, b, MaxCount)
    else
        return "Inviável"
    end
end

##

## Fase 2 para Problema de Produção:

# Seja o problema de produção para:
    #max {CT*X| A*X<=b, X>=0}

linear_programming(β, η, A, C, b, MaxCount)

simplex_Programming(A,C,b,MaxCount)





## Testes passados!!!!!!!!!!!



#
# Id = diagm(0=>fill(1., size(A,1)))
# A = [A Id]
# Id_C = transpose(fill(0,size(Id,1)))
# C = [C Id_C]
#
# ## Inicializa o que será utilizado
# B = A[:,β]
# N = A[:,η]
# cB = C[:,β]
# cN = C[:,η]
# XB = B\b
# b_barra = B\b
#
#
# w = cB*inv(B)
# global C_reduzido = -w*N + cN
# global pass = 0
# for count in 1:size(C_reduzido,2)
#     if C_reduzido[count] .<= 0
#         global pass +=1
#     end
# end
# if pass == size(C_reduzido,2)
#     global pass = false
# else
#     global pass = true
# end
#
#
# ## Loop que verifica o Custo reduzido.
# ## Parte principal do Programa
#
# #while pass
# save = 0
# for counter in 1:size(C_reduzido,2)
#     temp = C_reduzido[counter]
#     if counter == 1
#         global j = counter
#         global save = temp
#     elseif temp > save
#         global j = counter
#         global save = temp
#     end
# end
#
# a = N[:,j]
# dB_linha = -B\a
# #θ = findmax(XB + θ*dB>=0)
# if all(dB_linha .>= zeros(size(dB_linha,2)) )
#     return ("ilimitado")
# else
#     global save = 0
#     for i in 1:size(XB,1)
#         temp = XB[i]/abs(dB_linha[i])
#         if i == 1
#             global k = 1
#             global save = temp
#         elseif  temp < save
#             global k = i
#         end # if
#     end #for
# end #if
#
# aux = η[j]
# η[j] = β[k]
# β[k] = aux
# B = A[:,β]
# XB = B\b
# N = A[:,η]
# cB = C[:,β]
# cN = C[:,η]
# C_reduzido = -w*N + cN
# global pass = 0
# for count in 1:size(C_reduzido,2)
#     if C_reduzido[count] .<= 0
#         global pass +=1
#     end
# end
# if pass == size(C_reduzido,2)
#     global pass = false
# else
#     global pass = true
# end
# #end # while
#
# println("Ótimo Encontrado!")
# Xβ= B\b_barra
# Xη=zeros(size(A,2)-size(Xβ,1))
# X=[Xβ;Xη]
# Z = C*X
