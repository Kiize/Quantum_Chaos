using Base.Threads
using LinearAlgebra
using GLMakie
using DelimitedFiles
using ProgressMeter
using SparseArrays
using ArnoldiMethod
include("modules/module_geometries.jl")
include("modules/helper.jl")


R = RectGeom(0.0, 1.0, 0.0, 2.0)
S = BunStadium(R)
R_out = rect_around(S)

# Grid

N = 70 

x = range(R_out.Lx_min, R_out.Lx_max, length=N)
y = range(R_out.Ly_min, R_out.Ly_max, length=N)

#= points_in_stad = [(xi, yi) for xi in x, yi in y if isin_stadium(S, xi, yi)]
points_outside = [(xi, yi) for xi in x, yi in y if !isin_stadium(S, xi, yi)]

fig = Figure(size=(1200, 800))
ax = Axis(fig[1,1])
scatter!(ax, points_in_stad, color=:red)
scatter!(ax, points_outside, color=:blue)
fig =#

# v_nm
#= 
V_matrix = zeros(N, N) # La matrice v_nm

println("Calcolo della matrice V_nm $(N) × $(N)...\n")

# Cicli nidificati per n e m
@showprogress @threads for n = 1:N   # problem with @distributed
    # Vettore delle autofunzioni phi_n sulla griglia (linearizzato, solo Region II)
    # linearizzato_phi_n = P[:, :, n][Mask .== 1] 
    # Si consiglia di ciclare e usare gli indici per chiarezza e per evitare grosse allocazioni.

    for m = n:N # Simmetria: calcola solo la triangolare superiore (H_nm = H_mn)
        v_nm = integration_on_II(R, S, n, m, x, y)

        V_matrix[n, m] = v_nm
        V_matrix[m, n] = v_nm # Simmetria
    end
end

println("\n Matrice V_nm calcolata.\n")

open("data/matrix_v_nm_N$(N).txt", "w") do io
    writedlm(io, V_matrix)
end
 =#

# Eigenvalues
#= 
R_out_billiard = RectBilliard(R_out.Lx_max - R_out.Lx_min, R_out.Ly_max - R_out.Ly_min, N, N) # Rectangular Billiard
ham = rect_laplacian(R_out_billiard) # Hamiltonian
k = 70  # Number of eigenvalues to compute

# Diagonalization using Arnoldi method.

decomp, history = partialschur(ham, nev=k, tol=1e-6, which=:SR)  
history
E_num, Ψ_vecs = partialeigen(decomp)

# Save results on file to not redo them.

open("data/eigenvalues_rect_N$(N).txt", "w") do io
    writedlm(io, E_num)
end =#

# Diagonalization

E_all = vec(readdlm("data/eigenvalues_rect_N$(N).txt"))
V_matrix = readdlm("data/matrix_v_nm_N$(N).txt")
V_0 = 1000.0
k = 70

H = Diagonal(E_all) + V_0 * V_matrix

decomp, history = partialschur(H, nev=k, tol=1e-6, which=:SR)  
history
E_num_S, Ψ_vecs_S = partialeigen(decomp)

#= open("data/eigenvalues_Stadium_N$(N).txt", "w") do io
    writedlm(io, E_num_S)
end

open("data/eigenfun_Stadium_N$(N).txt", "w") do io
    writedlm(io, Ψ_vecs_S)
end =#

# Plot
fig = Figure(size=(1200, 800))
ax = Axis(fig[1,1])
scatter!(ax, E_num_S, marker=:circle, label="Stadium Eigenvalues")