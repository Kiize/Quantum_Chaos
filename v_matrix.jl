using Base.Threads
using LinearAlgebra
using DelimitedFiles
using ProgressMeter
include("modules/module_geometries.jl")
include("modules/helper.jl")

# Input data.
R = RectGeom(0.0, 1.0, 0.0, 2.0)
S = BunStadium(R)
R_out = rect_around(S)
N = 8
M = N^2
V_matrix = zeros(M, M) # La matrice v_nm


# Grid
function v_matrix(V_matrix::Matrix{Float64}, N::Int)
    M = N^2
    x = range(R_out.Lx_min, R_out.Lx_max, length=N)
    y = range(R_out.Ly_min, R_out.Ly_max, length=N)

    # Sorted energies.

    raw_states = [(
        (nx/(R_out.Lx_max - R_out.Lx_min))^2 + (ny/(R_out.Ly_max - R_out.Ly_min))^2, 
        nx, 
        ny
    ) for nx = 1:N, ny = 1:N]

    flattened_states = sort(reshape(raw_states, :), alg=PartialQuickSort(M))

    # v_nm


    println("Calcolo della matrice V_nm $(M) × $(M)...\n")

    # Cicli nidificati per n e m
    @showprogress @threads for n = 1:M   # problem with @distributed
        # Vettore delle autofunzioni phi_n sulla griglia (linearizzato, solo Region II)
        # linearizzato_phi_n = P[:, :, n][Mask .== 1] 
        # Si consiglia di ciclare e usare gli indici per chiarezza e per evitare grosse allocazioni.

        for m = n:M # Simmetria: calcola solo la triangolare superiore (H_nm = H_mn)
            v_nm = integration_on_II(R, S, n, m, x, y, flattened_states)  # @profview, @time

            V_matrix[n, m] = v_nm
            V_matrix[m, n] = v_nm # Simmetria
        end
    end

    println("\n Matrice V_nm calcolata.\n")
    
    open("data/data_stadium/matrix_v_nm_M$(M).txt", "w") do io
        writedlm(io, V_matrix)
    end

end

v_matrix(V_matrix, N)