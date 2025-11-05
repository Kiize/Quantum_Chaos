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