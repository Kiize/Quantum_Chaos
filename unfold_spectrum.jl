using LinearAlgebra
using GLMakie
using DelimitedFiles
using ProgressMeter
using LsqFit
using Statistics
#include("modules/module_geometries.jl")
#include("modules/helper.jl")

using QuantumChaos
# Geometry.

R = RectGeom(0.0, 1.0, 0.0, 2.0)
S = BunStadium(R)
R_out = rect_around(S)

# Application of Weyl's law.

N = 4900
fig = Figure(size=(1200, 800))

# LsqFit for parameters Stadium.

@. model(x, p) = p[1] / (4 * pi) * x - p[2] * (4 * pi) * sqrt(x) + p[3]

N_S = N ÷ 10 * 8   # eigenvalues for stadium
E_num_S = vec(readdlm("data/data_stadium/eigenvalues_Stadium_k$(N_S).txt"))
Sdata = 1:1:length(E_num_S)
p0 = [1., 1., 1.]

fitS = curve_fit(model, E_num_S, Sdata, p0)
AS, LS, CS = coef(fitS)

# LsqFit for parameters Rectangle.

N_R = N # eigenvalues for rect
E_num_R = vec(readdlm("data/data_stadium/eigenvalues_rect_k$(N_R).txt"))
Rdata = 1:1:length(E_num_R)

fitR = curve_fit(model, E_num_R, Rdata, p0)
AR, LR, CR = coef(fitR)


function weyl_energy(S::BunStadium, E::Vector{Float64}, A, L)
    #A = area_stadium(S)
    #L = perimeter_stadium(S)

    N = A / (4 * pi) .* E - L / (4 * pi) .* sqrt.(E)

    return N
end

weyl_energy(R::RectGeom, E::Vector{Float64}, A, L) = begin
    #A = area_rect(R)
    #L = perimeter_rect(R)

    N = A / (4 * pi) .* E - L / (4 * pi) .* sqrt.(E)

    return N
end

function weyl_law(N::Int, fig::Figure, AS, LS, AR, LR)
    N_R = N # eigenvalues for rect
    N_S = N ÷ 10 * 8   # eigenvalues for stadium
    E_num_S = vec(readdlm("data/data_stadium/eigenvalues_Stadium_k$(N_S).txt"))
    E_num_R = vec(readdlm("data/data_stadium/eigenvalues_rect_k$(N_R).txt"))


    ϵ = weyl_energy(S, E_num_S, AS, LS)
    s = diff(ϵ)
    println("variance s = $(std(s))")

    ϵ2 = weyl_energy(R, E_num_R, AR, LR)
    s2 = diff(ϵ2)
    println("variance s2 = $(std(s2))")

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

weyl_law(N, fig, AS, LS, AR, LR)
#save("figs/hist_ene_stadium_N$(N)_fit.png", fig)



