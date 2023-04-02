using CSV,LinearAlgebra,GLPK,Plots


A = [2 1;1 2]
C = [4 3]
b = [4;4]

## Matrizes das Restrições
R = [-1 -1]
bR = [-1]
##

    ## Função principal que organiza a ida para as fases I e II
function simplex_Programming(A , C, b, R, bR)

    ## Análise::
    V_zero = zeros(size(b))
    n = size(A,2)
    m = size(A,1)


    ## Verifica se b >=0, dessa forma, sendo direcionado para fase II
    if (b >= V_zero)
        η_init = zeros(Int8,n)
        β_init = zeros(Int8,m)

        for i in 1:n
            η_init[i] = Int(i)
        end
        for i in 1:m
            β_init[i] = Int(n+i)
        end
        ## União das Matrízes e Preparação para as Fases::

        V_zero = zeros(size(R,2))

        if R != transpose(V_zero)
            ## Tamanho das Linhas e colunas de A
            i = size(A,1)
            j = size(A,2)
            ##

            ## Número de Colunas de b (=1)
            jb = size(b,2)

            ## Tamanho das colunas de R tem que ser igual a de A
            jR = size(R,2)
            ## Tamanho das colunas do b das restrições(bR) tem que ser igual a de b
            jbR = size(bR,2)
            ##

            ## Verifica se há algum erro de tamanho de matrízes
            if( j != jR || jb != jbR)
                A = [A ; R]
                b = [b ; bR]
                return "ERROR! Matrix Mismatch"
            end
        end

        return simplex_FaseII(β_init, η_init, A, C, b)
    end

    ## Caso vá para a FaseI
    return simplex_FaseI(A,C,b)
end


## Início Fase I
function simplex_FaseI(A,C,b)

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

    Z_faseI, β_faseI, η_faseI = simplex_FaseII( β_init, η_init, Ae, -w, b )
    tol = 10^-5
    # Em algum momento isso aqui deve ser feito:
    if (Z_faseI == 0 || Z_faseI > -tol)
        return simplex_FaseII( β_faseI, η_faseI, A, C, b)
    else
        return "Inviável"
    end
end

## Início Fase II
function simplex_FaseII(β_init , η_init , A , C , b)

    ## Inicializa o que será utilizado
    Id = diagm(0=>fill(1., size(A,1)))
    Id_C = transpose(fill(0,size(Id,1)))
    C = [C Id_C]
    CT = transpose(C)
    A = [A Id]


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
        y =inv(B)*N[:,k]
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
        β[l] = k

        MatrizB = [inv(B) y]
        sizeMB = size(MatrizB,1)

        ultima_col = MatrizB[:,size(MatrizB,2)]
        escolhido = ultima_col[l]

        ## realiza as mudanças na matriz [B^-1 y]
        if escolhido != 0
             escolhido = MatrizB[l,:]/escolhido
        end

        for d in 1:sizeMB
            ## sim aqui dava pra ter otimizado melhor o código
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
        cB = C[:,β]
        cN = C[:,η]

        ## contador para o número de vezes que o loop acontece
        new_count+=1

        ## Etapa final:
        w = cB*inv(B)
        C_reduzido = w*N - cN
        β_init = β
        η_init = η
    end # while

    Xβ= B\b_barra
    Xη=zeros(size(A,2)-size(Xβ,1))
    X=[Xβ;Xη]
    Z = C*X

    return Z, η_init, β_init

end
##
	#Teste
Z, β, η = simplex_Programming(A, C, b, R, bR)


## Letra B
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



## Letra C
