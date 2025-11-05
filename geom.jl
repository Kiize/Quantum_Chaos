using LinearAlgebra
using GLMakie
include("module_geometries.jl")

R = RectGeom(0.0, 1.0, 0.0, 2.0)
S = BunStadium(R)
R_out = rect_around(S)

x = R_out.Lx_min:0.001:R_out.Lx_max
y = R_out.Ly_min:0.001:R_out.Ly_max

points_in_stad = [(xi, yi) for xi in x, yi in y if isin_stadium(S, xi, yi)]
points_outside = [(xi, yi) for xi in x, yi in y if !isin_stadium(S, xi, yi)]

fig = Figure(size=(1200, 800))
ax = Axis(fig[1,1])
scatter!(ax, points_in_stad, color=:red)
scatter!(ax, points_outside, color=:blue)
fig