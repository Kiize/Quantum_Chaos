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

# Application of Weyl's law.

N = 4900
fig = Figure(size=(1200, 800))

function weyl_energy(S::BunStadium, E::Vector{Float64})
    A = area_stadium(S)
    L = perimeter_stadium(S)

    # TO DO: fit parameters

    N = A / (4 * pi) .* E - L / (4 * pi) .* sqrt.(E)

    return N
end

weyl_energy(R::RectGeom, E::Vector{Float64}) = begin
    A = area_rect(R)
    L = perimeter_rect(R)

    # TO DO: fit parameters

    N = A / (4 * pi) .* E - L / (4 * pi) .* sqrt.(E)

    return N
end

function weyl_law(N::Int, fig::Figure)
    N_R = N # eigenvalues for rect
    N_S = N ÷ 10 * 8   # eigenvalues for stadium
    E_num_S = vec(readdlm("data/data_stadium/eigenvalues_Stadium_k$(N_S).txt"))
    E_num_R = vec(readdlm("data/data_stadium/eigenvalues_rect_k$(N_R).txt"))


    ϵ = weyl_energy(S, E_num_S)
    s = diff(ϵ)

    ϵ2 = weyl_energy(R, E_num_R)
    s2 = diff(ϵ2)
    ax = Axis(fig[1,1], xlabel = L"Spacing $s$", ylabel=L"P(s)", title="Distribution of the spacings in the Bunimovich stadium for k = $(N_S)")
    ax2 = Axis(fig[1,2], xlabel = L"Spacing $s$", ylabel=L"P(s)", title="Distribution of the spacings in the Rectangular billiard for k = $(N_R)")

    hist!(ax, s, normalization=:pdf)
    hist!(ax2, s2, normalization=:pdf)

    # PDFs.
    p(s) = exp(-s)
    goe(s) = pi * s / 2 * exp(-pi * s^2 / 4)

    s_plot = minimum(s):0.01:maximum(s)
    p_values = p.(s_plot)
    goe_values = goe.(s_plot)

    s2_plot = minimum(s2):0.01:maximum(s2)
    p2_values = p.(s2_plot)
    goe2_values = goe.(s2_plot)

    lines!(ax, s_plot, p_values, color = :red, label="Poisson")
    lines!(ax, s_plot, goe_values, color = :green, label="GOE")

    lines!(ax2, s2_plot, p2_values, color = :red, label="Poisson")
    lines!(ax2, s2_plot, goe2_values, color = :green, label="GOE")

    axislegend(ax)
    axislegend(ax2)
end

weyl_law(N, fig)
save("figs/hist_ene_stadium_N$(N).png", fig)
