using Arpack
using LinearAlgebra
using CSV
using DataFrames

H = CSV.File("Hamiltonian_100x100.csv", header=0) |> DataFrame |> Matrix
#@time X = eigvals(H)
@time X = eigs(H, nev=9999, ritzvec=false);
CSV.write("solved_100x100_2.csv",  Tables.table(X), writeheader=false)
