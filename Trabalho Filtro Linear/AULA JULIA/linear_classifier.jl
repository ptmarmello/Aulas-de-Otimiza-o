using DataFrames, CSV, JuMP, Plots, GLPK

##Set PATH, your current folder path..
PATH = pwd()
data_tumores = CSV.read(joinpath(PATH, "data_tumors.csv"))

tipo_tumor = data_tumores[:,2]
# tipo_tumor = data_tumores.Classification

caracteristicas = data_tumores[:,3:end]
carac_treino = caracteristicas[1:400,:]
carac_teste = caracteristicas[401:end,:]

μ = [8,10]

scatter(carac_treino[:,μ[1]],carac_treino[:,μ[2]])

function return_benigno_maligno_vector(caracteristicas, tipo_tumor)
    lin,col = size(caracteristicas)
    mal = Array{Any}(undef, 0,col)
    ben = Array{Any}(undef, 0,col)
    
    for i=1:length(caracteristicas[:,1])
        if tipo_tumor[i] == "M"
            mal = vcat(mal,Array(caracteristicas[i,:])')
        else
            ben = vcat(ben,Array(caracteristicas[i,:])')
        end
    end
    
    return mal, ben
end


mal, ben = return_benigno_maligno_vector(caracteristicas,tipo_tumor)


p = plot()
scatter(mal[:,μ[1]],mal[:,μ[2]])
scatter!(ben[:,μ[1]],ben[:,μ[2]])


model = Model(with_optimizer(GLPK.Optimizer))

@variable(model, a[μ])
@variable(model, c)
@variable(model, ϵ[i = 1:length(carac_treino[:,1])] >= 0)

for i in 1:length(carac_treino[:,1])
    if tipo_tumor[i] == "M"
        @constraint(model, sum(carac_treino[i, j]*a[j] for j in μ) + c >= -ϵ[i])
    elseif tipo_tumor[i] == "B"
        @constraint(model, sum(carac_treino[i, j]*a[j] for j in μ) + c <= ϵ[i]-1)
    end
end

@objective(model, Min, sum(ϵ))

optimize!(model)
termination_status(model)

a = value.(a)
c = value.(c)
obj_val = objective_value(model)

f(x) = -(a[μ[1]]*x)/a[μ[2]] - c/a[μ[2]]

p = plot()
scatter!(mal[:,μ[1]],mal[:,μ[2]])
scatter!(ben[:,μ[1]],ben[:,μ[2]])
plot!(f)

function get_vals(a)
    [a[i] for i in μ]
end

get_vals(a)