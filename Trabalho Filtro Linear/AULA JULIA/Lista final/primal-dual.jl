using JuMP, GLPK, LinearAlgebra,Plots, MathOptInterface

##1-] inicializa com X0>0 S0>0 p0 e k=0 viáveis

##2-] Teste de otimalidade: se S[k]' * x[k] < ϵ -> ótimo
# Passo caso não seja ótimo:


##3-] ρ[k] ∈(0,1] e definir:
    ## μ = ρ[k]*(S[k])' * x[k]/n
    ## X[k] = diagnm(x1[k],x2[k],x3[k]....xn[k])
    ## S[k] = diagnm(s1[k],s2[k],s3[k]....sn[k])
    #Resolver os sistemas lineares:
        ### solve: A*dX[k] =0
        ### solve: A'dy[k] + ds[k] =0
        ### solve: Sk*dx[k] + X[k]*ds[k]=μ[k]*e - X[k]*S[k]*e
##4-]
    # βP[k] = min{1,α min(-x[k][i]/dX[k][i])}, dX[k][i]<0
    # βD[k] = min{1,α min(-[Sk][i]/dS[k][i])}, dS[k][i]<0
#

#5-] encontrar:    X[k+1] = x[k] +βP[k]*dX[k]
#               p[k+1] = p[k] + βD[k]*dP[k]
#               s[k+1]= s[k] + βD[k]*dS[k]

#6-] acrescentar k:=k+1 e voltar ao passo 2.

#Função que inicia as matrizes que serão utilizadas para o Problema de Produção
function start()
    global A = [2 1;1 2]
    global C = [4 3]
    global b = [4;4]
    global p = 0
    global a = 0.999
    global MaxCount = 0
end

#Reinicia as matrizes e variáveis, caso seja necessario
function restart()
    start()
end

#Cria as matrizes iniciais dos problemas de forma "randomica"
function randStart(m,n)
    A = rand(Int8, m,n)
    C = rand(Int8, 1,n)
    b = (rand(Int8, m,1))
    p = 0.0001
    a = 0.999
    MaxCount = 0

    return A,C,b,p,a,MaxCount
end

## Função usada para chamar e usar o GLPK
function usingGLPK(A,C,b)
    (m,n) = size(A)
    modeloPrimal = Model(GLPK.Optimizer)
    set_optimizer_attribute(modeloPrimal, "msg_lev", 3)

    @variable(modeloPrimal, x[1:n])
    @constraint(modeloPrimal, (A*x)' .<= b)
    @constraint(modeloPrimal, x .>= 0)
    @objective(modeloPrimal, Max, sum(C*x))
    ansr, time_to_opt, gar1,gar2,gar3 = @timed begin
        optimize!(modeloPrimal)
    end
    time_by_solver = solve_time(modeloPrimal)
    status = termination_status(modeloPrimal)
    valor = value.(x)
    Z_Primal_GLPK = objective_value(modeloPrimal)
    has_duals(modeloPrimal)
    dual_status(modeloPrimal)
    Z_Dual_GLPK = dual_objective_value.(modeloPrimal)
    return valor,iter,Z_Primal_GLPK,Z_Dual_GLPK,time_to_opt,time_by_solver
end


## Funções para o Simplex
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


    w = zeros(n)
    w[n]= -1

    ### Call fase II para w

    Z_faseI, β_faseI, η_faseI = simplex_FaseII( β_init, η_init, Ae, -w', b, MaxCount)
    tol = 10^-5
    # Em algum momento isso aqui deve ser feito:
    if (Z_faseI == 0 || Z_faseI > -tol)
        return simplex_FaseII( β_faseI, η_faseI, A, C, b, MaxCount)
    else
        return "Inviável"
    end
end
function simplex_FaseII(β_init , η_init , A , C , b, MaxCount)

    Id = diagm(0=>fill(1., size(A,1)))
    A = [A Id]
    Id_C = zeros(size(C,2))
    C = [C Id_C']

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
    save_C_reduzido1 = zeros(1,1)
    save_C_reduzido2 = zeros(1,1)
    var = 0
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

        ##Falha de loop proposital

        # for i in 1:length(MaxCount)
        #     if C_reduzido == save_C_reduzido[i]
        #         break
        #     end
        # end
        # save_C_reduzido[MaxCount] = C_reduzido

        if MaxCount == 1
            save_C_reduzido1 = C_reduzido
        elseif MaxCount == 2
            aux = save_C_reduzido1
            save_C_reduzido1 = C_reduzido
            save_C_reduzido2 = aux
        elseif MaxCount >=3
            if C_reduzido == save_C_reduzido2
                Xβ= B\b_barra
                Xη=zeros(size(A,2)-size(Xβ,1))
                X=[Xβ;Xη]
                Z = C*X
                var +=1
                if var > 1
                    return(MaxCount, X, Z,"Ilimitado")
                end
            end
        end

    end # while

    Xβ= B\b_barra
    Xη=zeros(size(A,2)-size(Xβ,1))
    X=[Xβ;Xη]
    Z = C*X

    return MaxCount,X, Z, [β;η]
end

## Primeira Questão:::: Criar a função Pontos Interiores
## Função dos Pontos Interiores (resolve o Primal Dual)
function primalDual(A,b,C,p,a, MaxCount)


 ## Etapa 1
    (m,n) = size(A)
    Indx = 1:n
    Indy = n+1:n+m
    Inds = n+m+1:n+m+n

    #Inicialização
    x = ones(n) #x>=0
    y = ones(m)
    s = ones(n) #s>=0

    t = ones(n)
    # z = [x;y;s]
    μ = p*x'*s/n
 ## Etapa 2
    gap = x'*s
    Ε = ones(length(x),1)

    global Tx=1.
    global Ts=1.

    while all(s'*x > 0.000001)
        Fk = [A'*y+s-C';A*x-b;Diagonal(x)*Diagonal(s)*Ε - Ε*μ]
        Jk = [zeros(n,n) A' diagm(0=>fill(1., n));A zeros(m,m) zeros(m,n); Diagonal(s) zeros(n,m) Diagonal(x)]

        d = -Jk\Fk
        dx = d[Indx]
        dy = d[Indy]
        ds = d[Inds]

        #Passo implementado



        global saveTx = 0
        global saveTs = 0
        for i = 1:n
            if d[Indx[i]] <0
                global saveTx = -x[i]/dx[i]
                if i==1
                    global Tx = saveTx
                elseif saveTx < Tx
                    Tx = saveTx
                end #if
            end #if
        end #for

        for j in 1:n
            if (ds[j] <0)
                global saveTs = -s[j]/ds[j]
                if j==1
                    global Ts = saveTs
                elseif saveTs < Ts
                    global Ts = saveTs
                end
            end #if
        end # for
        Tx = min(1,Tx*a)
        Ts = min(1,Ts*a)

        x = x + Tx*dx
        y = y + Ts*dy
        s = s + Ts*ds

        # z = [x;y;s]
        μ = p*x'*s/n
        # gap = [gap x'*s]
        MaxCount = MaxCount +1

        ## Verificar aqui a convergência, infinito, ilimitado, inviável
        t_temp = b-A*x
        u_temp = C'-A'*y-s
        if t_temp == 0 || u_temp == 0
            if dx>0 && C*dx<0
                return ("Ilimitado!")
            elseif ds>0 && b*dy>0
                return ("Ilimitado!")
            end # if
        end
        if MaxCount > 2000
            return ("Ilimitado")
        end
    end # while
    Z_Primal = C*x
    Z_Dual = y'*b



    return MaxCount, Z_Primal,Z_Dual,x
end


start()
@timed begin
    primalDual(A,b,C,p,a, MaxCount)
end # begin
@elapsed primalDual(A,b,C,p,a, MaxCount)

## Letra b


##Teste para o problema de Produção usando GLPK
restart()
usingGLPK(A,C,b)

##Teste Para o problema de produção usando Simplex

resultado_Simplex, time_simplex, lixo1, lixo2,lixo3 = @timed begin
    simplex_Programming(A,C,b,MaxCount)
end

##Teste Para o problema de produção usando Primal-Dual (Pontos Interiores)

resultado_PD, time_PD, lixo1, lixo2,lixo3 = @timed begin
    primalDual(A,b,C,p,a, MaxCount)
end


Z_Primal_GLPK
Z_Dual_GLPK
Z_Primal_PD

Z_Dual_PD
Z_Simplex


count_Simplex
count_PD

time_GLPK = zeros(3)
X_GLPK = zeros(3)
counter_simplex = zeros(3)
time_Simplex = zeros(3)
X_Simplex = zeros(3)
counter_PD = zeros(3)
time_pd = zeros(3)
X_pd = zeros(3)

for i in 1: 3
    m=2
    n=2
    A,C,b,p,a,MaxCount = randStart(m,n)
    resposta_GLPK = @timed begin
        usingGLPK(A,C,b)
    end

    params, time_GLPK[i], lixo1,lixo2,lixo3 = resposta_GLPK
    X_GLPK, z_primal_glpk,z_dual_glpk, time1,time2 = params


    resposta_Simplex = @timed begin
        simplex_Programming(A,C,b,MaxCount)
    end

    params, time_Simplex[i],lixo1,lixo2,lixo3 = resposta_Simplex
    counter_simplex[i], X_Simplex, outros = params


    resposta_PD = @timed begin
        primalDual(A,b,C,p,a, MaxCount)
    end

    params, time_pd[i],lixo1,lixo2,lixo3 = resposta_PD
    counter_PD[i], Z_primal_pd,Z_dal_pd,X_pd = params
end
plot(time_pd, counter_PD)
plot(time_Simplex, counter_simplex)
