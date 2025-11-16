using Base.Threads
using LinearAlgebra
using GLMakie
using DelimitedFiles
using ProgressMeter
using SparseArrays
using ArnoldiMethod
include("modules/module_geometries.jl")
include("modules/helper.jl")

R_vec = readdlm("data/data_stadium/rect_around.txt")   # Need to convert the vector to a RectGeom.
R_out = RectGeom(R_vec...)

# Eigenvalues of the outer Rectangular billiard R_out.
# N is the number of discretized points, while k is the number of eigenvalues we want to evaluate using the Arnoldi method.
function num_eigen_rect(R_out::RectGeom, k::Int)
    N = Int(sqrt(k))
    R_out_billiard = RectBilliard(R_out.Lx_max - R_out.Lx_min, R_out.Ly_max - R_out.Ly_min, N, N) # Rectangular Billiard
    ham = rect_laplacian(R_out_billiard) # Hamiltonian N^2 × N^2

    println("\n Diagonalizing with k=$(k)...\n")
    # Diagonalization using Arnoldi method.
    decomp, history = partialschur(ham, nev=k, tol=1e-6, which=:SR)  
    history
    E_num, _ = partialeigen(decomp)

    # Save results on file to not redo them.
    println("\n Saving results on file...\n")

    open("data/data_stadium/eigenvalues_rect_k$(k).txt", "w") do io
        writedlm(io, E_num)
    end
end

k_list = [1600, 2500, 3600, 4900]
#k_list = [64]
for k in k_list
    num_eigen_rect(R_out, k)
end