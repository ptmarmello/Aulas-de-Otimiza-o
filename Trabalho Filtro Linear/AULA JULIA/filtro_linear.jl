#using Pkg
#Pkg.add("CSV")
#Pkg.add("JuMP")
#Pkg.add("GLPK")
#Pkg.add("DataFrames")
using CSV, JuMP, GLPK, DataFrames, Plots
tabela_ons = CSV.read("D:\\Users\\user\\Desktop\\otimizacao\\serie_carga_ons.csv")
carga = tabela_ons[:,4]
## Questão 1
## Definir um modelo de otimização:
model = Model(GLPK.Optimizer)
## Variável para o problema:
T = size(carga)[1]
@variable(model, τ[i = 1:T])
@variable(model, δ[i =1:T])
@variable(model, Δ[i =2:T-1])

for i in 1:T
    @constraint(model, δ[i] >= carga[i] - τ[i])
    @constraint(model, δ[i] >= -(carga[i]) + τ[i])
end
for i in 2:(T-1)
    @constraint(model, Δ[i] >= (τ[i+1] - τ[i]) - (τ[i] - τ[i-1]))
    @constraint(model, Δ[i] >= -(τ[i+1] - τ[i]) + (τ[i] - τ[i-1]))
end
λ = 1
obj = sum(δ) + (λ * (sum(Δ)))
@objective(model, Min, obj)
optimize!(model)
status = termination_status(model)
vτ = value.(τ)

## Questão 2
x = 1:T
plot(x,carga)

λ = [1,10,100,1000,10000]
obj1 = sum(δ) + (λ[1] * (sum(Δ)))
obj2 = sum(δ) + (λ[2] * (sum(Δ)))
obj3 = sum(δ) + (λ[3] * (sum(Δ)))
obj4 = sum(δ) + (λ[4] * (sum(Δ)))
obj5 = sum(δ) + (λ[5] * (sum(Δ)))


@objective(model, Min, obj1)
optimize!(model)
status = termination_status(model)
vτ1 = value.(τ)

@objective(model, Min, obj2)
optimize!(model)
status = termination_status(model)
vτ2 = value.(τ)

@objective(model, Min, obj3)
optimize!(model)
status = termination_status(model)
vτ3 = value.(τ)

@objective(model, Min, obj4)
optimize!(model)
status = termination_status(model)
vτ4 = value.(τ)

@objective(model, Min, obj5)
optimize!(model)
status = termination_status(model)
vτ5 = value.(τ)

plot(x,vτ1)
plot(x,vτ2)
plot(x,vτ3)
plot(x,vτ4)
plot(x,vτ5)


modelhp = Model(GLPK.Optimizer)
@variable(modelhp, τ[i = 1:T])
@variable(modelhp, t[i = 2:T-1])
for i in 2:T-1
    @constraint(modelhp, τ[i] = t[1])
end
obj = sum((carga-τ)^2) + (λ * (sum()))
@objective(modelhp, Min, obj)
optimize!(modelhp)
status = termination_status(modelhp)
vτ = value.(τ)

## Questão 3
