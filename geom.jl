using Base.Threads
using LinearAlgebra
using GLMakie
using DelimitedFiles
using ProgressMeter
using SparseArrays
using ArnoldiMethod
#include("modules/module_geometries.jl")
#include("modules/helper.jl")

using QuantumChaos


R = RectGeom(0.5, 1.0, 0.5, 2.0)
S = BunStadium(R)
R_out = rect_around(S)

# Write R_out on file to export the same configuration.
#= 
open("data/data_stadium/rect_around.txt", "w") do io
    writedlm(io, RectGeom_as_vec(R_out))
end
 =#
# Grid

N = 200 
k = N  # Number of eigenvalues to compute

x = range(R_out.Lx_min, R_out.Lx_max, length=N)
y = range(R_out.Ly_min, R_out.Ly_max, length=N)

# Plot stadium.
points_in_stad = [(xi, yi) for xi in x, yi in y if isin_stadium(S, xi, yi)]
points_outside = [(xi, yi) for xi in x, yi in y if !isin_stadium(S, xi, yi)]
points_circ = [(xi, yi) for xi in x, yi in y if isin_circle(S.circ_left, xi, yi) || isin_circle(S.circ_right, xi, yi)]

fig = Figure(size=(1200, 800))
ax = Axis(fig[1,1], xlabel="x", ylabel="y")
scatter!(ax, points_in_stad, color=:red)
scatter!(ax, points_outside, color=:blue)
#scatter!(ax, points_circ, color=:green)
fig

# Plot.
#= 
fig = Figure(size=(1200, 800))
ax = Axis(fig[1,1], xlabel = "Energy index", ylabel="relative error")#, title="Relative error of Stadium Eigenvalues")
ax2 = Axis(fig[1,2], xlabel = "Energy index", ylabel="E")#, title="Relative error of Stadium Eigenvalues")


function rel_err_plot(N::Int, ax::Axis)
    E_num_S = vec(readdlm("data/data_stadium/eigenvalues_Stadium_k$(N).txt"))
    rel_E = relative_error(E_num_S)
    
    # Plot
    scatter!(ax, rel_E, label="N = $(N)") 
end

function stadium_E_plot(N::Int, ax::Axis)
    E_num_S = vec(readdlm("data/data_stadium/eigenvalues_Stadium_k$(N).txt"))
    
    # Plot
    scatter!(ax, E_num_S, label="N = $(N)") 
end

N_list = [1600, 2500, 3600, 4900]
#N_list = [50]
for n in N_list
    n = n ÷ 2
    rel_err_plot(n, ax)
    stadium_E_plot(n, ax2)
end

axislegend(ax)
axislegend(ax2)

save("figs/stadium_eigenvalues_compare.png", fig) =#
