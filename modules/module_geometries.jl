using LinearAlgebra

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

rect_eigenfun(R, k::Int, x, y) = begin
    base = base_map(R, k) 
    (_, nx, ny) = base[k]
    return rect_eigenfun(R, nx, ny, x, y)
end

#= function eigenfun_on_mesh(R::RectGeom, k::Int, x_vals::StepRangeLen{Float64}, y_vals::StepRangeLen{Float64})
    psi_vals = [rect_eigenfun(R, k, x, y) for x in x_vals, y in y_vals]
    return psi_vals
end =#

# Returns v_{nm}       
function integration_on_II(R::RectGeom, S::BunStadium, n::Int, m::Int, x_vals::StepRangeLen{Float64}, y_vals::StepRangeLen{Float64})
    dx = step(x_vals)
    dy = step(y_vals)

    integral = 0.0

    for x in x_vals
        for y in y_vals
            if !isin_stadium(S, x, y)
                psi_n = rect_eigenfun(R, n, x, y)
                psi_m = rect_eigenfun(R, m, x, y)
                integral += psi_n * psi_m * dx * dy
            end
        end
    end

    return integral
end