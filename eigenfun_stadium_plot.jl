using LinearAlgebra
using GLMakie
using DelimitedFiles
using ProgressMeter
include("modules/module_geometries.jl")
include("modules/helper.jl")

# Geometry:
# S is built from the rectangle R and R_out is the surrounding rectangle.

R = RectGeom(0.0, 1.0, 0.0, 2.0)
S = BunStadium(R)
R_out = rect_around(S)

# Parameters:
# N is the number of points along each axis, k is the number of eigenvalues we evaluated previously.

N = 70
M = N^2
k = M ÷ 10 * 8

# Sorted states flattened_states = ϕₘ.

raw_states = [(
    (nx/(R_out.Lx_max - R_out.Lx_min))^2 + (ny/(R_out.Ly_max - R_out.Ly_min))^2, 
    nx, 
    ny
) for nx = 1:N, ny = 1:N]

flattened_states = sort(reshape(raw_states, :), alg=PartialQuickSort(M))

# Coefficients cₘ.
coeff = readdlm("data/data_stadium/eigenfun_Stadium_k$(k).txt") 

# Figures.
f = Figure()
f2 = Figure()

# Heatmap wavefunction. We want to plot the ik-th eigenstates.
ik = 1
plot_stadium_eigenstate(k, R_out, N, f, f2, flattened_states, coeff, ik)