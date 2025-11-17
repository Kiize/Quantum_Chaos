using Base.Threads
using LinearAlgebra
using GLMakie
using DelimitedFiles
using ProgressMeter
include("modules/module_geometries.jl")
include("modules/helper.jl")


R = RectGeom(0.0, 1.0, 0.0, 2.0)
S = BunStadium(R)
R_out = rect_around(S)

N = 40
M = N^2
k = M ÷ 10 * 8

# Sorted energies.

raw_states = [(
    (nx/(R_out.Lx_max - R_out.Lx_min))^2 + (ny/(R_out.Ly_max - R_out.Ly_min))^2, 
    nx, 
    ny
) for nx = 1:N, ny = 1:N]

flattened_states = sort(reshape(raw_states, :), alg=PartialQuickSort(M))

# P matrix.

function calc_P_matrix(BaseMap, x_grid, y_grid, Lx, Ly)
    M = length(BaseMap)
    Nx, Ny = length(x_grid), length(y_grid)
    N_points = Nx * Ny
    P = zeros(Float64, N_points, M) # P è N_points x M
    A = sqrt(4.0 / (Lx * Ly)) # Fattore di normalizzazione

    @showprogress "Calculating P matrix..." for m = 1:M
        # Recupera i numeri quantici (nx, ny)
        E_m, nx_m, ny_m = BaseMap[m] 
        
        # Pre-calcolo Vettoriale 1D (Veloce)
        Sx = sin.(nx_m * pi * x_grid / Lx) 
        Sy = sin.(ny_m * pi * y_grid / Ly)

        # Prodotto Esterno (Nx x Ny) e appiattimento in P[:, m]
        Z_m = A * Sx * Sy'
        P[:, m] .= vec(Z_m) 
    end

    # Implementa la Maschera (Mask) per forzare psi=0 fuori dal biliardo
    # (Devi implementare una funzione per questo. Assumo che 'is_in_stadium' esista)
    Mask_linear = zeros(N_points)
    for (k, (x, y)) in enumerate(zip(repeat(x_grid, outer=Nx), repeat(y_grid, inner=Ny)))
        if isin_stadium(S, x, y) # Funzione esterna, da definire!
            Mask_linear[k] = 1.0
        end
    end
    return P, Mask_linear
end




# Heatmap wavefunction.


coeff = readdlm("data/data_stadium/eigenfun_Stadium_k$(k).txt")
#x = range(R_out.Lx_min, R_out.Lx_max, length=N)
#y = range(R_out.Ly_min, R_out.Ly_max, length=N)
#R_out_billiard = RectBilliard(R_out.Lx_max - R_out.Lx_min, R_out.Ly_max - R_out.Ly_min, N, N) # Rectangular Billiard
f = Figure()
f2 = Figure()

function eigenfun_old(xvec, yvec, flattened_states::Vector{Tuple{Float64, Int64, Int64}}, coeff::Matrix{Float64}, k::Int, R_out) 
    coeff2D = reshape(coeff, N, N, k)
    ψ = zeros(N, N)
    for i in 1:k
        for j in eachindex(xvec)
            for h in eachindex(yvec)
                ψ[j, h] = coeff2D[j, h, i] * rect_eigenfun(R_out, xvec[j], yvec[h], flattened_states[i])
            end
        end
    end
    
    return ψ
end

function eigenfun(coeff, ik::Int, x_grid, y_grid, Lx, Ly)
    c_k = coeff[:, ik]
    N_grid = length(x_grid)
    # Calcola P una sola volta!
    P, Mask_linear = calc_P_matrix(flattened_states, x_grid, y_grid, Lx, Ly)

    # Moltiplicazione Matrice-Vettore: P * c_k
    # Psi_k_linear è il vettore N_points x 1
    Psi_k_linear = P * c_k 
    
    # Applicazione della Maschera e rimodellamento
    Psi_k_masked = Psi_k_linear #.* Mask_linear
    return reshape(Psi_k_masked, N_grid, N_grid)
end
# Function to plot the k-th eigenstate, it also saves the plots in the "figs" folder.
function plot_stadium_eigenstate(k, R::RectGeom, N::Int, f::Figure, f2::Figure, flattened_states, coeff)
    Lx = R.Lx_max - R.Lx_min
    Ly = R.Ly_max - R.Ly_min
    
    hx = Lx / (N + 1)
    hy = Ly / (N + 1) 
    
    # Points for the heatmap and contour plot.
    x_int = range(R.Lx_min + hx, R.Lx_max - hx, length=N)
    y_int = range(R.Ly_min + hy, R.Ly_max - hy, length=N)
    
    # Reshape it in 2D array for plotting.
    #psi_2D = eigenfun_old(x_int, y_int, flattened_states, coeff, k, R)
    psi_2D = eigenfun(coeff, 3, x_int, y_int, Lx, Ly)
    #psi_2D = reshape(psi_vector, N, N)

    # Heatmap.
    ax = Axis(f[1, 1], title="Probability density |ψ($k)|^2 (FDM)", xlabel="x", ylabel="y")
    hm = heatmap!(ax, x_int, y_int, abs2.(psi_2D))
    Colorbar(f[1, 2], hm)
    # Contour.
    ax2 = Axis(f2[1, 1], title="Probability density |ψ($k)|^2 (FDM)", xlabel="x", ylabel="y")
    cmap =  :diverging_bkr_55_10_c35_n256
    contourf!(ax2, x_int, y_int, abs2.(psi_2D), colormap=cmap)
    #contour3d!(ax2, x_int, y_int, log.(abs2.(psi_2D)); levels = 14, colormap = cmap, transparency = true, linewidth = 5) # 3d contour on the right.

    # Save figures.
    save("figs/figs_stadium/rect_eigenstate_$k.png", f)
    save("figs/figs_stadium/rect_eigenstate_$(k)_cont.png", f2)
    #f, f2
    display(f)
    #display(f2)
end

plot_stadium_eigenstate(k, R_out, N, f, f2, flattened_states, coeff)