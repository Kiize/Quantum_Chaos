using LinearAlgebra
using SparseArrays
using ArnoldiMethod
using Kronecker
using GLMakie
using DelimitedFiles

using QuantumChaos


# Definition of the Rectangular Billiard B, the Laplacian operator ham over B and the number of eigenvalues k to compute. 

B = RectBilliard(1.0, sqrt(2), 200, 300) # Rectangular Billiard
ham = rect_laplacian(B) # Hamiltonian
k = 100  # Number of eigenvalues to compute

# Write B on file to export the same configuration.

open("data/rect_billiard.txt", "w") do io
    writedlm(io, RectBilliard_as_vec(B))
end

# Diagonalization using Arnoldi method.

decomp, history = partialschur(ham, nev=k, tol=1e-6, which=:SR)  
history
E_num, Ψ_vecs = partialeigen(decomp)

# Save results on file to not redo them.

open("data/eigenvalues.txt", "w") do io
    writedlm(io, E_num)
end

open("data/eigenvecs.txt", "w") do io
    writedlm(io, Ψ_vecs)
end

# Plot k-th eigenstate and compare numerical and analytical eigenvalues.

plot_eigenstate(k, Ψ_vecs, B)
compare_energies(E_num, B, k)