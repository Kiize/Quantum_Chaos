using LinearAlgebra
using DelimitedFiles
#using ArnoldiMethod
using Arpack
using ProgressMeter
include("modules/module_geometries.jl")
include("modules/helper.jl")

# Diagonalization of the M × M Hamiltonian matrix
# H = E_all 𝕀_M + V_0 V_matrix.
# We take the first k eigenvalues.
function stadium_diag(M::Int, k::Int)
    E_all = vec(readdlm("data/data_stadium/eigenvalues_rect_k$(M).txt"))
    V_matrix = readdlm("data/data_stadium/matrix_v_nm_M$(M).txt")
    V_0 = 10 * E_all[end]

    H = Diagonal(E_all) + V_0 * V_matrix

    println("\n Diagonalizing with M = $(M), k = $(k)...\n")

    #= decomp, history = partialschur(H, nev=k, tol=1e-6, which=:SR)  
    history
    E_num_S, Ψ_vecs_S = partialeigen(decomp) =#

    E_num_S, Ψ_vecs_S = eigs(H, nev=k, which=:SM)

    println("\n Saving results on file...\n")

    open("data/data_stadium/eigenvalues_Stadium_k$(k).txt", "w") do io
        writedlm(io, E_num_S)
    end

    open("data/data_stadium/eigenfun_Stadium_k$(k).txt", "w") do io
        writedlm(io, Ψ_vecs_S)
    end
end

#M_list = [10, 20, 50, 70, 100]
M_list = [1600, 2500, 3600, 4900]
for M in M_list
    k = M ÷ 2
    stadium_diag(M, k)
end