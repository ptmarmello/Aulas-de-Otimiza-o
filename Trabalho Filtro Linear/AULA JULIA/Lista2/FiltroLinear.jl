#using Pkg
# Pkg.add("JuMP")
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("Plots")
# Pkg.add("GLPK")
#Pkg.add("Ipopt")

using DataFrames, CSV, JuMP, Plots, GLPK, Ipopt
PATH=pwd()
#-----------------------------------------------------------------------------#
#Questão 1:
#-----------------------------------------------------------------------------#
# Formulando o que foi apresentado no PDF

modelo = Model(GLPK.Optimizer)
T=size(carga)[1]
#Variáveis:
@variable(modelo,τ[i =1:T])
@variable(modelo, δ[i =1:T])
@variable(modelo, Δ[i =2:T-1])

for i in 1:T
    @constraint(modelo,δ[i] >= carga[i] - τ[i])
    @constraint(modelo,δ[i] >= -(carga[i]) + (τ[i]))
end
for i in 2:(T-1)
    @constraint(modelo,Δ[i] >= (τ[i+1] - τ[i]) - (τ[i]-τ[i-1]) )
    @constraint(modelo,Δ[i] >= -(τ[i+1] - τ[i]) + (τ[i]-τ[i-1]) )
end

#-----------------------------------------------------------------------------#
#Questão 2:
#-----------------------------------------------------------------------------#
# Leitura do arquivo ::
tabela_cargaOns = CSV.read("serie_carga_ons.csv")
y = tabela_cargaOns[:,4]
T=size(y)[1]
x=1:T
plot(x,y)


#Variáveis:
modelo = Model(GLPK.Optimizer)
@variable(modelo,τ[i =1:T])
@variable(modelo, δ[i =1:T])
@variable(modelo, Δ[i =2:T-1])

for i in 1:T
    @constraint(modelo,δ[i] >= y[i] - τ[i])
    @constraint(modelo,δ[i] >= -(y[i]) + (τ[i]))
end
for i in 2:(T-1)
    @constraint(modelo,Δ[i] >= (τ[i+1] - τ[i]) - (τ[i]-τ[i-1]) )
    @constraint(modelo,Δ[i] >= -(τ[i+1] - τ[i]) + (τ[i]-τ[i-1]) )
end


λ = 1
obj = sum(δ) + (λ* (sum(Δ)))
@objective(modelo, Min, obj)
optimize!(modelo)
status = termination_status(modelo)

##Plot:
X=value.(τ)
plot(X,label="λ=1")
## Aplicando o modelo com diferentes λ:::

λ = 10
obj1 = sum(δ) + (λ* (sum(Δ)))
@objective(modelo, Min, obj1)
optimize!(modelo)
status = termination_status(modelo)
X1=value.(τ)
plot(X1,label="λ=10")

λ = 200
obj2 = sum(δ) + (λ* (sum(Δ)))
@objective(modelo, Min, obj2)
optimize!(modelo)
status = termination_status(modelo)
X2=value.(τ)
plot(X2,label="λ=200")

λ = 1200
obj3 = sum(δ) + (λ* (sum(Δ)))
@objective(modelo, Min, obj3)
optimize!(modelo)
status = termination_status(modelo)
X3=value.(τ)
plot(X3,label="λ=1200")

λ = 10000
obj4 = sum(δ) + (λ* (sum(Δ)))
@objective(modelo, Min, obj4)
optimize!(modelo)
status = termination_status(modelo)
X4=value.(τ)
plot(X4,label="λ=10000")

## Modelo HP

modeloHP = Model(Ipopt.Optimizer)
@variable(modeloHP,τ[i =1:T])

λ = 1
objHP = sum((y[i] - τ[i])^2 for i in 1:T) + λ* sum(((τ[i+1]-τ[i]) - (τ[i]-τ[i-1]))^2 for i in 2:T-1)
@objective(modeloHP, Min, objHP)
optimize!(modeloHP)
status = termination_status(modeloHP)
X=value.(τ)
plot(X)

λ = 20
objHP2 = sum((y[i] - τ[i])^2 for i in 1:T) + λ* sum(((τ[i+1]-τ[i]) - (τ[i]-τ[i-1]))^2 for i in 2:T-1)
@objective(modeloHP, Min, objHP2)
optimize!(modeloHP)
status = termination_status(modeloHP)
X=value.(τ)
plot(X)

λ = 150
objHP3 = sum((y[i] - τ[i])^2 for i in 1:T) + λ* sum(((τ[i+1]-τ[i]) - (τ[i]-τ[i-1]))^2 for i in 2:T-1)
@objective(modeloHP, Min, objHP3)
optimize!(modeloHP)
status = termination_status(modeloHP)
X=value.(τ)
plot(X)

λ = 2000
objHP4 = sum((y[i] - τ[i])^2 for i in 1:T) + λ* sum(((τ[i+1]-τ[i]) - (τ[i]-τ[i-1]))^2 for i in 2:T-1)
@objective(modeloHP, Min, objHP4)
optimize!(modeloHP)
status = termination_status(modeloHP)
X=value.(τ)
plot(X)

λ = 10000
objHP5 = sum((y[i] - τ[i])^2 for i in 1:T) + λ* sum(((τ[i+1]-τ[i]) - (τ[i]-τ[i-1]))^2 for i in 2:T-1)
@objective(modeloHP, Min, objHP5)
optimize!(modeloHP)
status = termination_status(modeloHP)
X=value.(τ)
plot(X)

#-----------------------------------------------------------------------------#
#Questão 3:
#-------------------------------------------------------------------------------
tabela = CSV.read("FinancialDistress.csv")
y = tabela[:,13]
T=size(y)[1]
x = 1:T
plot(x,y, label = "original")

cargaComRuido = zeros(T)
for i in 1:T
    cargaComRuido[i] = y[i] + rand(-0.5:0.1:0.5)
end

plot!(x,cargaComRuido,label = "Ruido")

#Variáveis:
NewModelo = Model(GLPK.Optimizer)
@variable(NewModelo,τ[i =1:T])
@variable(NewModelo, δ[i =1:T])
@variable(NewModelo, Δ[i =2:T-1])

@variable(NewModelo,ε[i =1:T])
@variable(NewModelo, ϕ[i =1:T])
@variable(NewModelo, Φ[i =2:T-1])


## para a original
for i in 1:T
    @constraint(NewModelo,δ[i] >= y[i] - τ[i])
    @constraint(NewModelo,δ[i] >= -(y[i]) + (τ[i]))
end
for i in 2:(T-1)
    @constraint(NewModelo,Δ[i] >= (τ[i+1] - τ[i]) - (τ[i]-τ[i-1]) )
    @constraint(NewModelo,Δ[i] >= -(τ[i+1] - τ[i]) + (τ[i]-τ[i-1]) )
end

##Para o ruído
for i in 1:T
    @constraint(NewModelo,ϕ[i] >= cargaComRuido[i] - ε[i])
    @constraint(NewModelo,ϕ[i] >= -(cargaComRuido[i]) + (ε[i]))
end
for i in 2:(T-1)
    @constraint(NewModelo,Φ[i] >= (ε[i+1] - ε[i]) - (ε[i]-ε[i-1]) )
    @constraint(NewModelo,Φ[i] >= -(ε[i+1] - ε[i]) + (ε[i]-ε[i-1]) )
end

λ = 100
obj = sum(δ) + (λ* (sum(Δ)))
@objective(NewModelo, Min, obj)
optimize!(NewModelo)
status = termination_status(NewModelo)
X=value.(τ)
plot(X,label = "original 100")



objRand = sum(ϕ) + (λ* (sum(Φ)))
@objective(NewModelo, Min, objRand)
optimize!(NewModelo)
status = termination_status(NewModelo)
X=value.(ε)
plot!(X,label="random 100")



## Discrepância
cargaComRuido = zeros(T)
for i in 1:T
    cargaComRuido[i] = y[i] + rand(-5:0.1:5)
end
for i in 1:T
    @constraint(NewModelo,ϕ[i] >= cargaComRuido[i] - ε[i])
    @constraint(NewModelo,ϕ[i] >= -(cargaComRuido[i]) + (ε[i]))
end
for i in 2:(T-1)
    @constraint(NewModelo,Φ[i] >= (ε[i+1] - ε[i]) - (ε[i]-ε[i-1]) )
    @constraint(NewModelo,Φ[i] >= -(ε[i+1] - ε[i]) + (ε[i]-ε[i-1]) )
end

λ = 100
obj = sum(δ) + (λ* (sum(Δ)))
@objective(NewModelo, Min, obj)
optimize!(NewModelo)
status = termination_status(NewModelo)
X=value.(τ)
plot(X,label = "original 100")

objRand = sum(ϕ) + (λ* (sum(Φ)))
@objective(NewModelo, Min, objRand)
optimize!(NewModelo)
status = termination_status(NewModelo)
X=value.(ε)
plot!(X,label="random 100")
