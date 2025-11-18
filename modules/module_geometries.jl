using LinearAlgebra, Integrals

struct RectGeom
    Lx_min::Float64
    Lx_max::Float64
    Ly_min::Float64
    Ly_max::Float64
end

struct CircGeom
    radius::Float64
    x_center::Float64
    y_center::Float64
end

struct BunStadium 
    rect::RectGeom
    circ_right::CircGeom
    circ_left::CircGeom
end

function BunStadium(R::RectGeom)
    radius = (R.Ly_max + R.Ly_min)/2.0
    circ_right = CircGeom(radius, R.Lx_max, radius)
    circ_left = CircGeom(radius, R.Lx_min, radius)
    return BunStadium(R, circ_right, circ_left)
end

function area_stadium(S::BunStadium)
    area_rect = (S.rect.Lx_max - S.rect.Lx_min) * (S.rect.Ly_max - S.rect.Ly_min)
    area_circ = S.circ_left.radius^2 * pi

    return area_rect + area_circ
end

function area_rect(R::RectGeom)
    return (R.Lx_max - R.Lx_min) * (R.Ly_max - R.Ly_min)
end

function perimeter_stadium(S::BunStadium)
    per_circ = 2 * pi * S.circ_left.radius
    per_rect = 2 * (S.rect.Lx_max - S.rect.Lx_min)

    return per_circ + per_rect
end

function perimeter_rect(R::RectGeom)
    return 2 * (R.Lx_max - R.Lx_min) + 2 * (R.Ly_max - R.Ly_min)
end

function isin_rect(R::RectGeom, x::Float64, y::Float64)
    return (R.Lx_min <= x <= R.Lx_max) && (R.Ly_min <= y <= R.Ly_max)
end

function isin_circle(C::CircGeom, x::Float64, y::Float64)
    return (x - C.x_center)^2 + (y - C.y_center)^2 <= C.radius^2
end

function isin_stadium(S::BunStadium, x::Float64, y::Float64)
    in_rect = isin_rect(S.rect, x, y)
    in_right_circle = isin_circle(S.circ_right, x, y)
    in_left_circle = isin_circle(S.circ_left, x, y)
    return in_rect || in_right_circle || in_left_circle
end

function rect_around(minx::Float64, maxx::Float64, miny::Float64, maxy::Float64)
    return RectGeom(minx, maxx, miny, maxy)
end

rect_around(S::BunStadium) = rect_around(S.rect.Lx_min - S.circ_left.radius,
                                            S.rect.Lx_max + S.circ_right.radius,
                                            S.rect.Ly_min,
                                            S.rect.Ly_max)

# Function to export the RectGeom parameters as a vector, used to write it on file.
function RectGeom_as_vec(R::RectGeom)
    r = Union{Float64, Float64, Float64, Float64}[R.Lx_min::Float64, R.Lx_max::Float64, R.Ly_min::Float64, R.Ly_max::Float64]
    return r
end

# Eigenfunctions of the external rectangle.
function rect_eigenfun(R::RectGeom, nx::Int, ny::Int, x::Float64, y::Float64)
    Lx = R.Lx_max - R.Lx_min
    Ly = R.Ly_max - R.Ly_min

    psi_x = sqrt(2/Lx) * sin(pi/Lx * nx * x)
    psi_y = sqrt(2/Ly) * sin(pi/Ly * ny * y)

    return psi_x * psi_y    
end

# base_map maps the pair (nx, ny) to one single index k. 
# N is the maximum quantum number considered in each direction and also such that the Hamiltonian matrix will be of size N × N.
function base_map(R::RectGeom, N::Int)
    # Evaluate the energy and use it to sort the states.
    Lx = R.Lx_max - R.Lx_min
    Ly = R.Ly_max - R.Ly_min

    raw_states = []

    for nx = 1:N
        for ny = 1:N
            E = (nx/Lx)^2 + (ny/Ly)^2
            push!(raw_states, (E, nx, ny))
        end
    end

    # 2. Ordina e Tronca a M stati
    #M = 300 # La dimensione della matrice H sarà 300x300
    sort!(raw_states, by = x -> x[1])

    # BaseMap è la lista finale M x (E, nx, ny)
    return raw_states[1:N]
end

rect_eigenfun(R, x, y, flattened_states_sorted) = begin
    # assume flattened_states_sorted is already the k-th element 
    #base = flattened_states_sorted[1:k]
    base = flattened_states_sorted
    #base = base_map(R, k) 
    #(_, nx, ny) = base[k]
    (_, nx, ny) = base
    return rect_eigenfun(R, nx, ny, x, y)
end

#= function eigenfun_on_mesh(R::RectGeom, k::Int, x_vals::StepRangeLen{Float64}, y_vals::StepRangeLen{Float64})
    psi_vals = [rect_eigenfun(R, k, x, y) for x in x_vals, y in y_vals]
    return psi_vals
end =#

# Returns v_{nm}       
function integration_on_II(R::RectGeom, S::BunStadium, n::Int, m::Int, x_vals::StepRangeLen{Float64}, y_vals::StepRangeLen{Float64}, flattened_states::Vector{Tuple{Float64, Int64, Int64}})
    # Old code, faster.
    dx = step(x_vals)
    dy = step(y_vals)

    integral = 0.0

    for x in x_vals
        for y in y_vals
            if !isin_stadium(S, x, y)
                psi_n = rect_eigenfun(R, x, y, flattened_states[n])
                psi_m = rect_eigenfun(R, x, y, flattened_states[m])
                integral += psi_n * psi_m * dx * dy
            end
        end
    end

    # New code using Integrals.
    #= 
    f(u, p) = rect_eigenfun(R, n, u[1], u[2]) * rect_eigenfun(R, m, u[1], u[2]) * !isin_stadium(S, u[1], u[2])
    domain = ([x_vals[1], y_vals[1]], [x_vals[end], y_vals[end]])
    prob = IntegralProblem(f, domain)
    sol = solve(prob, HCubatureJL(); reltol = 1e-3, abstol = 1e-3)

    integral = sol.u =#

    return integral
end

function relative_error(E::Vector{Float64})
    E_out = similar(E)
    E_out[1] = 1
    for i in 2:length(E)
        E_out[i] = abs((E[i] - E[i - 1]))/E[i - 1]
    end

    return E_out
end

# P matrix. The P matrix contains the eigenfunctions of the rectangular billiard evaluated on the grid points.
# Given the billiard of length Lx × Ly we can discretize it and construct the grid x_grid × y_grid. 
function calc_P_matrix(flattened_states, x_grid, y_grid, Lx, Ly)
    len = length(flattened_states)
    Nx, Ny = length(x_grid), length(y_grid) # In principle we could implement the case where Nx ≠ Ny.
    M = Nx * Ny
    P = zeros(Float64, M, len) # P is M × len.
    A = sqrt(4.0 / (Lx * Ly)) # Normalization factor.

    @showprogress "Calculating P matrix..." for m = 1:len
        _, nx_m, ny_m = flattened_states[m] 
        
        # Eigenfunctions on x and y.
        Sx = sin.(nx_m * pi * x_grid / Lx) 
        Sy = sin.(ny_m * pi * y_grid / Ly)

        Z_m = A * Sx * Sy'
        P[:, m] .= vec(Z_m) 
    end

    # Mask to force ψ (P_matrix) = 0 outside the stadium.
    Mask_linear = zeros(M)
    for (k, (x, y)) in enumerate(zip(repeat(x_grid, outer=Nx), repeat(y_grid, inner=Ny)))
        if isin_stadium(S, x, y) 
            Mask_linear[k] = 1.0
        end
    end
    return P, Mask_linear
end

# Function to evaluate the state ψ = ∑ₘ cₘ ϕₘ, where cₘ = coeff are the eigenstates evaluated in stadium_diag.jl.
function eigenfun(coeff, flattened_states, ik::Int, x_grid, y_grid, Lx, Ly)
    c_k = coeff[:, ik]
    N_grid = length(x_grid)
    P, Mask_linear = calc_P_matrix(flattened_states, x_grid, y_grid, Lx, Ly)

    # Psi_k_linear is a vector of length M.
    Psi_k_linear = P * c_k 
    
    # Masking, if needed, our vector.
    Psi_k_masked = Psi_k_linear #.* Mask_linear
    return reshape(Psi_k_masked, N_grid, N_grid)
end


# Function to plot the ik-th eigenstate, it also saves the plots in the "figs" folder.
function plot_stadium_eigenstate(k, R::RectGeom, N::Int, f::Figure, f2::Figure, flattened_states, coeff, ik::Int)
    Lx = R.Lx_max - R.Lx_min
    Ly = R.Ly_max - R.Ly_min
    
    hx = Lx / (N + 1)
    hy = Ly / (N + 1) 
    
    # Points for the heatmap and contour plot.
    x_int = range(R.Lx_min + hx, R.Lx_max - hx, length=N)
    y_int = range(R.Ly_min + hy, R.Ly_max - hy, length=N)
    
    # Reshape it in 2D array for plotting.
    #psi_2D = eigenfun_old(x_int, y_int, flattened_states, coeff, k, R)
    psi_2D = eigenfun(coeff, flattened_states, ik, x_int, y_int, Lx, Ly)
    #psi_2D = reshape(psi_vector, N, N)

    # Heatmap.
    ax = Axis(f[1, 1], title="Probability density |ψ($(ik))|^2 (FDM)", xlabel="x", ylabel="y")
    hm = heatmap!(ax, x_int, y_int, abs2.(psi_2D))
    Colorbar(f[1, 2], hm)
    # Contour.
    #ax2 = Axis(f2[1, 1], title="Probability density |ψ($(ik))|^2 (FDM)", xlabel="x", ylabel="y")
    ax2 = Axis3(f2[1,1]; aspect = (1,1,0.7), perspectiveness = 0.5)

    cmap = :viridis # :diverging_bkr_55_10_c35_n256
    contourf!(ax2, x_int, y_int, abs2.(psi_2D), colormap=cmap)
    contour3d!(ax2, x_int, y_int, abs2.(psi_2D); levels = 14, colormap = cmap, transparency = true, linewidth = 5) # 3d contour on the right.

    # Save figures.
    save("figs/figs_stadium/rect_eigenstate_$k.png", f)
    save("figs/figs_stadium/rect_eigenstate_$(k)_cont.png", f2)
    #f, f2
    display(f)
    #display(f2)
end