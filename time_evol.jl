using DelimitedFiles
using LinearAlgebra
using GLMakie
include("help_evol.jl")
#include("helper.jl")


# Example usage (replace indices as necessario)
E_all = vec(readdlm("eigenvalues.txt"))
Ψ_all = readdlm("eigenvecs.txt")    # matrix with eigenvectors as columns
B_vec = readdlm("rect_billiard.txt")
B = RectBilliard(B_vec...)

number_states = 5
c = fill(1/√number_states, number_states)                   # real coefficients (will be cast to Complex)
t = 0.0:0.1:1.0

# take first 3 eigenstates for the superposition
pe = PsiEvol(c, E_all[1:number_states], Ψ_all[:, 1:number_states], t)

# psi at time index 5
ψ_t5 = psi_at_index(pe, 5)

# psi at arbitrary time value
ψ_0_37 = psi_at_time(pe, 0.37)

# Time evolution plot
T = 10.0
fig = Figure(size = (1200,800))
ax = Axis(fig[1, 1], title="Evoluzione Temporale della Densità di Probabilità", xlabel="x", ylabel="y")
ax2 = Axis3(fig[1,2]; aspect = (1,1,0.7), perspectiveness = 0.5)
axs = [ax, ax2]
y = psi_at_time(pe, 0.0) 

hm = plot_at_t!(y, B, axs)

framerate = 20
timestamps = range(0, 1, step=1/framerate)

record(fig, "time_animation.mp4", 0.0:0.1:T; framerate = framerate) do t
    y = psi_at_time(pe, t) 
    hm = plot_at_t!(y, B, axs)
end