# Stage3_Vortex_Excitations/src/loaders.jl
#
# Stage3 needs no optimization of its own — R_match and R0 are loaded
# from Stage1 and Stage2 respectively, since the vortex term doesn't
# change the AA/AB physics those parameters were optimized for.

function load_Rmatch(stage1_dir::String, N_half::Int, L::Float64)::Float64
    path = result_path(stage1_dir, "Rmatch_optimum", (N=N_half, L=L))
    if !isfile(path)
        error("No saved R_match for N=$N_half, L=$L. Run Stage1's main_VMC.jl first.")
    end
    return load_run(path; run=1).result.R_opt
end

function load_R0(stage2_dir::String, num_part::Int, nr0_sq::Float64, h::Float64)::Float64
    path = result_path(stage2_dir, "R0_optimum", (N=num_part, nr0sq=nr0_sq, h=h))
    if !isfile(path)
        error("No saved R0 for N=$num_part, nr0sq=$nr0_sq, h=$h. Run Stage2's main_VMC.jl first.")
    end
    return load_run(path; run=1).result.R0_opt
end