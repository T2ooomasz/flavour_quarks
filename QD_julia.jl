using LinearAlgebra
using CSV
using DataFrames

site_size = 3000
site_elements = 130
PBC = "False"
V_max = 1 
r_0 = 1000
m_eff = 0.067
title =   "_" * string(site_elements) * "_" * string(site_size) * "_" * string(r_0) * "_" *  string(V_max) * "_" * string(PBC) * "_" * string(m_eff) * ".csv"

H = CSV.File("Hamiltonian" * title, header=0) |> DataFrame |> Matrix
println("starting calculation")
@time F = eigen(H)
CSV.write("eigenvalues_julia" * title,  Tables.table(F.values), writeheader=false)
CSV.write("psi_julia" * title,  Tables.table(F.vectors), writeheader=false)
