using LinearAlgebra
using GLMakie
include("helper.jl")

struct PsiEvol
    c::Vector{ComplexF64}        # coefficients of eigenstates
    E::Vector{Float64}          # eigenenergies
    psi_vecs::Matrix{ComplexF64}# eigenvectors (columns = eigenstates)
    t::Vector{Float64}          # time grid
    psi_t::Matrix{ComplexF64}   # psi at each time (columns = times)

    function PsiEvol(c_in, E_in, psi_in, t_in)
        c = ComplexF64.(vec(c_in))
        E = Float64.(vec(E_in))
        psi_vecs = ComplexF64.(psi_in)
        t = Float64.(collect(t_in))

        Nspace, Nk = size(psi_vecs)
        if length(c) != Nk || length(E) != Nk
            throw(ArgumentError("Mismatch: length(c)=$(length(c)), length(E)=$(length(E)), number of eigenvectors=$(Nk)"))
        end

        Nt = length(t)
        psi_t = zeros(ComplexF64, Nspace, Nt)

        # Efficient assembly: psi(t_j) = psi_vecs * (c .* exp(-im * E * t_j))
        for j in 1:Nt
            phases = c .* exp.(-im .* E .* t[j])
            psi_t[:, j] = psi_vecs * phases
        end

        new(c, E, psi_vecs, t, psi_t)
    end
end

# accessors
psi_at_index(pe::PsiEvol, j::Integer) = pe.psi_t[:, j]
psi_at_time(pe::PsiEvol, tval::Real) = pe.psi_vecs * (pe.c .* exp.(-im .* pe.E .* float(tval)))

function plot_at_t!(psi_vector_t, B::RectBilliard, axs::Array{Makie.AbstractAxis, 1})
    # Rimappa il vettore 1D alla griglia 2D (Nx x Ny)
    # Usiamo la funzione di reshaping (nota che devi scegliere l'ordine corretto)
    psi_2D = reshape(psi_vector_t, B.Nx, B.Ny)

    # Definisci il range dei punti (senza i bordi)
    x_int = range(B.hx, B.Lx - B.hx, length=B.Nx)
    y_int = range(B.hy, B.Ly - B.hy, length=B.Ny)

    heatmap!(axs[1], x_int, y_int, abs2.(psi_2D))

    #contourf!(axs[2], x_int, y_int, abs2.(psi_2D); levels = 14, colormap = :viridis)
    contour3d!(axs[2], x_int, y_int, abs2.(psi_2D); levels = 14, colormap = :viridis,
        transparency = true, linewidth = 5)
end