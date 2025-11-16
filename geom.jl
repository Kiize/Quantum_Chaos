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

open("data/rect_around.txt", "w") do io
    writedlm(io, RectGeom_as_vec(R_out))
end

# Grid

N = 200 
k = N  # Number of eigenvalues to compute

x = range(R_out.Lx_min, R_out.Lx_max, length=N)
y = range(R_out.Ly_min, R_out.Ly_max, length=N)

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

# Application of Weyl's law.

N = 4900 ÷ 2
fig = Figure(size=(1200, 800))

function weyl_energy(S::BunStadium, E::Vector{Float64})
    A = area_stadium(S)
    L = perimeter_stadium(S)

    # TO DO: fit parameters

    N = A / (4 * pi) .* E - L / (4 * pi) .* sqrt.(E)

    return N
end

function weyl_law(N::Int, fig::Figure)
    E_num_S = vec(readdlm("data/data_stadium/eigenvalues_Stadium_k$(N).txt"))


    ϵ = weyl_energy(S, E_num_S)
    s = diff(ϵ)
    ax = Axis(fig[1,1], xlabel = L"Spacing $s$", ylabel=L"P(s)", title="Distribution of the spacings for k = $(N)")


    hist!(ax, s, normalization=:pdf)

    # PDFs.
    p(s) = exp(-s)
    goe(s) = pi * s / 2 * exp(-pi * s^2 / 4)

    s_plot = minimum(s):0.01:maximum(s)
    p_values = p.(s_plot)
    goe_values = goe.(s_plot)

    lines!(ax, s_plot, p_values, color = :red, label="Poisson")
    lines!(ax, s_plot, goe_values, color = :green, label="GOE")
end


weyl_law(N, fig)
save("figs/hist_ene_stadium_N$(N).png", fig)
