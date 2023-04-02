using Pkg
Pkg.add("JuMP")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Ipopt")
Pkg.add("GLPK")

using DataFrames, JuMP, CSV, Ipopt, GLPK

PATH = pwd()
