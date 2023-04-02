#using Pkg
#Pkg.add("GLPK")
#Pkg.add("JuMP")
using JuMP,GLPK


## QUESTÃO 3::::
## DISCLAIMER:: A questão 3 não foi finalizada, por isso está tudo comentado.


## QUESTÃO 2: Letra a

model=Model(GLPK.Optimizer)

@variable(model, q)
@variable(model, g1)
@variable(model, g2)
@variable(model, g3)
@variable(model, f12)
@variable(model, f13)
@variable(model, f23)
c1=100
c2=150
c3=200
G1=5
G2=20
G3=12
F12=20
F13=20
F23=5

@constraint(model, q>=0)
@constraint(model, q<=15)

@constraint(model, g1>=0)
@constraint(model, g1<=G1)
@constraint(model, g2>=0)
@constraint(model, g2<=G2)
@constraint(model, g3>=0)
@constraint(model, g3<=G3)

@constraint(model, f12<=F12)
@constraint(model, f12>=-F12)
@constraint(model, f13<=F13)
@constraint(model, f13>=-F13)
@constraint(model, f23<=F23)
@constraint(model, f23>=-F23)

# @constraint(model, g1+g2+g3>=15)
@constraint(model, g1-f12-f13>=0)
@constraint(model, g1-f12-f13<=0)
@constraint(model, g2+f12-f23>=0)
@constraint(model, g2+f12-f23<=0)
@constraint(model, g3+f12+f23>=0)
@constraint(model, g3+f12+f23<=0)

obj2 = c1*g1 +c2*g2+c3*g3 + (15-q)*5000

@objective(model, Min, obj2)
optimize!(model)
status = termination_status(model)
qvalue= value.(q)
g1value=value.(g1)
g2value=value.(g2)
g3value=value.(g3)

answer= c1*g1value+ c2*g2value+c3*g3value+(15-qvalue)*5000


##Questão 2: Letra b


model2=Model(GLPK.Optimizer)

@variable(model2, q)
@variable(model2, g1)
@variable(model2, g2)
@variable(model2, g3)
@variable(model2, f12)
@variable(model2, f13)
@variable(model2, f23)
c1=100
c2=150
c3=200
G1=5
G2=20
G3=12
F12=20
F13=20
F23=5

@constraint(model2, q>=0)
@constraint(model2, q<=15)

@constraint(model2, g1>=0)
@constraint(model2, g1<=G1)
@constraint(model2, g2>=0)
@constraint(model2, g2<=G2)
@constraint(model2, g3>=0)
@constraint(model2, g3<=G3)

@constraint(model2, f12<=F12)
@constraint(model2, f12>=-F12)
@constraint(model2, f13<=F13)
@constraint(model2, f13>=-F13)
@constraint(model2, f23<=F23)
@constraint(model2, f23>=-F23)

@constraint(model2, g1+g2+g3>=15)
# @constraint(model2, g1-f12-f13>=0)
# @constraint(model2, g1-f12-f13<=0)
# @constraint(model2, g2+f12-f23>=0)
# @constraint(model2, g2+f12-f23<=0)
# @constraint(model2, g3+f12+f23>=0)
# @constraint(model2, g3+f12+f23<=0)

obj3 = c1*g1 +c2*g2+c3*g3 + (15-q)*5000

@objective(model2, Min, obj3)
optimize!(model2)
status = termination_status(model)
qvalue= value.(q)
g1value=value.(g1)
g2value=value.(g2)
g3value=value.(g3)

answer2= c1*g1value+ c2*g2value+c3*g3value+(15-qvalue)*5000


##QUESTÃO 3::: letra e
# 
# modelo = Model(GLPK.Optimizer)
# @variable(modelo,N)
# n=10
# @variable(modelo, c[1:10])
# @variable(modelo, G[1:10])
# @variable(modelo, Rup[1:10])
# @variable(modelo, Rdown[1:10])
# @variable(modelo, d[1:24])
# @variable(modelo, g[1:10, 1:24])
# @variable(modelo, Δg[1:10, 1:24])
# @variable(modelo, fbat[1:24])
# @variable(modelo, ΔGbat[1:24])
# @variable(modelo, Gbat)
# @variable(modelo, bat[1:24])
#
#
# for i=1:10
#     c[i].==2*i
#     G[i].==22-2*i
#     Rup[i].==i
#     Rdown[i].==i
# end
# for t=1:24
#     d[t].==60*(1+sin(t/12))
# end
# C1=100
# C2=50
# Gbat=8
# X=1000
# # i=1:10
# # t=1:24
#
# @constraint(modelo, sum(sum(g) + fbat) <= d + ΔGbat)
# @constraint(modelo, g[i, t]>=0)
# @constraint(modelo, g[i, t]<=G[i])
# @constraint(modelo, Δg[i, t]=g[i][t+1]-g[i][t])
# @constraint(modelo, Δg[i, t]<=Rup[i])
# @constraint(modelo, Δg[i, t]>=-Rdown[i])
# @constraint(modelo, ΔGbat[t]=Gbat-bat[t])
# @constraint(modelo, bat[t]>=0)
# @constraint(modelo, bat[t]<=Gbat)
# @constraint(modelo, fbat[t]>=0)
# @constraint(modelo, fbat[t]<=bat[t])
#
# obj =  N*sum( sum(c[i]*g[i,t], i=1:10) + F(d[t]), t=1:24)
# @objective(model, Min, obj2)
# optimize!(model)
# status = termination_status(model)
#
# ans=value.(N)
