using LinearAlgebra
using GLMakie
include("helper.jl")    # Import the RectBilliard struct.

# Define the evolution of a generic state psi_t.
# We write a generic state psi_t as a superposition of the eigenstates psi_vecs: psi_t = ∑ᵢ c[i] psi_vecs[i].
# We then now that eigenvectors evolve in time with a phase given by their energy E.
struct PsiEvol
    c::Vector{ComplexF64}        # coefficients of eigenstates.
    E::Vector{Float64}          # eigenenergies.
    psi_vecs::Matrix{ComplexF64}# eigenvectors (columns = eigenstates).
    t::Vector{Float64}          # time grid.
    psi_t::Matrix{ComplexF64}   # psi at each time (columns = times).

    # Constructor to evaluate psi_t given the coefficients, the energies, the eigenvectors and the time.
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

# Accessors to find psi_t at a specific time index or a specific time.
psi_at_index(pe::PsiEvol, j::Integer) = pe.psi_t[:, j]
psi_at_time(pe::PsiEvol, tval::Real) = pe.psi_vecs * (pe.c .* exp.(-im .* pe.E .* float(tval)))

# Function to get the x, y coordinates of the billiard for the heatmap.
function shape_rect(B::RectBilliard)

    # Definisci il range dei punti (senza i bordi)
    x_int = range(B.hx, B.Lx - B.hx, length=B.Nx)
    y_int = range(B.hy, B.Ly - B.hy, length=B.Ny)

    return x_int, y_int
end