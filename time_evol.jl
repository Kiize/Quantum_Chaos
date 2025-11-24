using DelimitedFiles
using LinearAlgebra
using GLMakie
using ProgressMeter
#include("modules/help_evol.jl")


# Read the stored eigenenergies, eigenstates and the rectangular billiard.
E_all = vec(readdlm("data/eigenvalues.txt"))
Ψ_all = readdlm("data/eigenvecs.txt")    
B_vec = readdlm("data/rect_billiard.txt")   # Need to convert the vector to a RectBilliard.
B = RectBilliard(B_vec...)

# Construct the initial state.
number_states = 7   # the first eigenstates considered.
c = fill(1/√number_states, number_states)   # equal probaility for each eigenstates.
T = 1.0    # final time.
time = 0.0:0.001:T

# Our state.
pe = PsiEvol(c, E_all[1:number_states], Ψ_all[:, 1:number_states], time)

# Time evolution plot.
fig = Figure(size = (1200,800))
ax = Axis(fig[1, 1], title="Time evolution of probaility density", xlabel="x", ylabel="y")
ax2 = Axis3(fig[1,2]; aspect = (1,1,0.7), perspectiveness = 0.5)
axs = [ax, ax2]

y = psi_at_time(pe, 0.0) # initial state.
x_int, y_int = shape_rect(B)    # x, y for heatmap plot.
psi_2D = reshape(y, B.Nx, B.Ny) # we want to plot the amplitude |ψ|^2.

hm = heatmap!(axs[1], x_int, y_int, abs2.(psi_2D))  # heatmap on the left.
ct = contourf!(axs[2], x_int, y_int, abs2.(psi_2D); levels = 14, colormap = :viridis)   # basis of the 3d contour on the right.
ct3d = contour3d!(axs[2], x_int, y_int, abs2.(psi_2D); levels = 14, colormap = :viridis,
        transparency = true, linewidth = 5) # 3d contour on the right.


framerate = 60
# Record animation.
record(fig, "figs/time_animation.mp4"; framerate = framerate) do io # we want to use showprogress to see how much time does it take to record.
    @showprogress for t in time
        y = psi_at_time(pe, t) 
        psi_2D = reshape(y, B.Nx, B.Ny)

        # Update the plots.
        hm[3] = abs2.(psi_2D) # update, not make new one.
        ct[3] = abs2.(psi_2D)
        ct3d[3] = abs2.(psi_2D)
        autolimits!.(axs)
        recordframe!(io)    # record the frame
    end
end