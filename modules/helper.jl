using LinearAlgebra, SparseArrays, Kronecker, GLMakie

# Definition of the Rectangular Billiard with dimension Lx × Ly and number of points Nx × Ny. 
struct RectBilliard 
    Lx::Float64
    Ly::Float64
    Nx::Int
    Ny::Int
    hx::Float64 
    hy::Float64 
end

# Automatic calculation of the steps hx and hy.
function RectBilliard(Lx::Float64, Ly::Float64, Nx::Int, Ny::Int)
    hx = Lx / (Nx + 1)
    hy = Ly / (Ny + 1)
    
    return RectBilliard(Lx, Ly, Nx, Ny, hx, hy)
end

# Constructor to convert floating values of Nx, Ny in Int.
RectBilliard(Lx::Float64, Ly::Float64, Nx::Float64, Ny::Float64) = RectBilliard(Lx, Ly, Int(Nx), Int(Ny))   # converte Nx, Ny in Int

# Function to export the RectBilliard parameters as a vector, used to write it on file.
function RectBilliard_as_vec(B::RectBilliard)
    b = Union{Float64, Float64, Int, Int}[B.Lx::Float64, B.Ly::Float64, B.Nx::Int, B.Ny::Int]
    return b
end

# Function to construct the Laplacian operator over the Rectangular Billiard B using Finite Difference Method (FDM).
function rect_laplacian(B::RectBilliard)
    # Laplacian matrix in 1D 
    Ax = spdiagm(
        0 => fill(-2.0, B.Nx), 
        -1 => fill(1.0, B.Nx - 1), 
        1 => fill(1.0, B.Nx - 1)
    )

    Ay = spdiagm(
        0 => fill(-2.0, B.Ny), 
        -1 => fill(1.0, B.Ny - 1), 
        1 => fill(1.0, B.Ny - 1)
    )

    # Identity matrix
    I_Nx = spdiagm(0 => fill(1.0, B.Nx)) 
    I_Ny = spdiagm(0 => fill(1.0, B.Ny)) 


    # To construct the 2D Laplacian Matrix we need to use Kronecker products (⊗, kron). We use kron to ensure the resulting matrix is sparse. 
    Hx = -(1 / (B.hx^2)) * kron(Ax, I_Ny)
    Hy = -(1 / (B.hy^2)) * kron(I_Nx, Ay)

    # H = -∇²
    H = Hx + Hy

    return H
end

# Function to plot the k-th eigenstate, it also saves the plots in the "figs" folder.
function plot_eigenstate(k, eigen_vecs, B::RectBilliard)
    println("--- Numerical result for (Nx=$(B.Nx), Ny=$(B.Ny)) ---")
    println("Numerical energy E_num[$k]: $(E_num[k])")

    # Select the k-th eigenvector.
    psi_vector = eigen_vecs[:, k]
    # Reshape it in 2D array for plotting.
    psi_2D = reshape(psi_vector, B.Nx, B.Ny)

    # Points for the heatmap and contour plot.
    x_int = range(B.hx, B.Lx - B.hx, length=B.Nx)
    y_int = range(B.hy, B.Ly - B.hy, length=B.Ny)

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
    save("figs/rect_eigenstate_$k.png", f)
    save("figs/rect_eigenstate_$(k)_cont.png", f2)
    f, f2
end

# Function to compare numerical and analytical eigenvalues, it also saves the comparison plot in the "figs" folder.
function compare_energies(energies, B::RectBilliard, k)
    energies_ana = []
    for nx in 1:k
        for ny in 1:k
            E_ana = pi^2 * ( (nx/B.Lx)^2 + (ny/B.Ly)^2 ) # Analytical formula for the rectangular billiard
            push!(energies_ana, E_ana)
        end
    end

    energies_ana = sort(energies_ana)[1:k]

    # Plot comparison.
    f = Figure()
    ax = Axis(f[1, 1], xlabel="E_num", ylabel="E_ana")
    scatter!(ax, energies, energies_ana, label="Energies")
    lines!(ax, energies_ana[1]:energies_ana[end], energies_ana[1]:energies_ana[end], color=:red, label="y=x", linewidth=1)
    axislegend(position = :rb)
    save("figs/rect_energy_comparison_$k.png", f)
    f
end