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

# Write R_out on file to export the same configuration.
#= 
open("data/data_stadium/rect_around.txt", "w") do io
    writedlm(io, RectGeom_as_vec(R_out))
end

# Grid

N = 200 
k = N  # Number of eigenvalues to compute

x = range(R_out.Lx_min, R_out.Lx_max, length=N)
y = range(R_out.Ly_min, R_out.Ly_max, length=N)
 =#
# Plot stadium.
#= points_in_stad = [(xi, yi) for xi in x, yi in y if isin_stadium(S, xi, yi)]
points_outside = [(xi, yi) for xi in x, yi in y if !isin_stadium(S, xi, yi)]

fig = Figure(size=(1200, 800))
ax = Axis(fig[1,1])
scatter!(ax, points_in_stad, color=:red)
scatter!(ax, points_outside, color=:blue)
fig =#

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

# Heatmap wavefunction.
N = 50
k = 2000
ψ = readdlm("data/data_stadium/eigenfun_Stadium_k$(k).txt")
#x = range(R_out.Lx_min, R_out.Lx_max, length=N)
#y = range(R_out.Ly_min, R_out.Ly_max, length=N)
#R_out_billiard = RectBilliard(R_out.Lx_max - R_out.Lx_min, R_out.Ly_max - R_out.Ly_min, N, N) # Rectangular Billiard

# Function to plot the k-th eigenstate, it also saves the plots in the "figs" folder.
function plot_stadium_eigenstate(k, eigen_vecs, R::RectGeom, N::Int)
    Lx = R.Lx_max - R.Lx_min
    Ly = R.Ly_max - R.Ly_min
    
    hx = Lx / (N + 1)
    hy = Ly / (N + 1) 
    # Select the k-th eigenvector.
    psi_vector = eigen_vecs[:, k]
    # Reshape it in 2D array for plotting.
    psi_2D = reshape(psi_vector, N, N)

    # Points for the heatmap and contour plot.
    x_int = range(R.Lx_min + hx, R.Lx_max - hx, length=N)
    y_int = range(R.Ly_min + hy, R.Ly_max - hy, length=N)

    # Heatmap.
    f = Figure()
    ax = Axis(f[1, 1], title="Probability density |ψ($k)|^2 (FDM)", xlabel="x", ylabel="y")
    heatmap!(ax, x_int, y_int, abs2.(psi_2D))

    # Contour.
    f2 = Figure()
    ax2 = Axis(f2[1, 1], title="Probability density |ψ($k)|^2 (FDM)", xlabel="x", ylabel="y")
    cmap =  :diverging_bkr_55_10_c35_n256
    contourf!(ax2, x_int, y_int, abs2.(psi_2D), colormap=cmap)

    # Save figures.
    save("figs/figs_stadium/rect_eigenstate_$k.png", f)
    save("figs/figs_stadium/rect_eigenstate_$(k)_cont.png", f2)
    #f, f2
    display(f)
    display(f2)
end

plot_stadium_eigenstate(1, ψ, R_out, N)