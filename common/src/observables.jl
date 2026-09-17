# common/src/observables.jl

"""
    accumulate_gr!(gr_histogram, x_coord, y_coord, L)

In-place accumulation of the intra-species pair correlation function g(r)
for a single MC configuration. `gr_histogram` should be preallocated once
outside the MC loop (zeros(Float64, num_bins)) and passed in; call once
per sampled (decorrelated) step. Call `normalize_gr!` once, after the loop,
to convert the raw histogram into g(r).
"""
function accumulate_gr!(gr_histogram::Vector{Float64}, x_coord::Vector{Float64},
                        y_coord::Vector{Float64}, L::Float64)
    num_part = length(x_coord)
    num_bins = length(gr_histogram)
    dr = (L/2) / num_bins

    @inbounds for i in 1:num_part
        for j in (i+1):num_part
            dx = get_periodic_difference(x_coord[i], x_coord[j], L)
            dy = get_periodic_difference(y_coord[i], y_coord[j], L)
            r = sqrt(dx^2 + dy^2)
            if r < L/2
                bin_index = Int(floor(r / dr)) + 1
                gr_histogram[bin_index] += 2
            end
        end
    end
    return nothing
end

"""
    accumulate_g_AB_r!(gr_histogram, x_coord_A, y_coord_A, x_coord_B, y_coord_B, L)

In-place accumulation of the inter-species (A-B) pair correlation function,
for two-component systems. Same call pattern as `accumulate_gr!`.
"""
function accumulate_g_AB_r!(gr_histogram::Vector{Float64}, x_coord_A::Vector{Float64},
                          y_coord_A::Vector{Float64}, x_coord_B::Vector{Float64},
                          y_coord_B::Vector{Float64}, L::Float64)
    num_part_A = length(x_coord_A)
    num_part_B = length(x_coord_B)
    num_bins = length(gr_histogram)
    dr = (L/2) / num_bins

    @inbounds for i in 1:num_part_A
        for j in 1:num_part_B
            dx = get_periodic_difference(x_coord_A[i], x_coord_B[j], L)
            dy = get_periodic_difference(y_coord_A[i], y_coord_B[j], L)
            r = sqrt(dx^2 + dy^2)
            if r < L/2
                bin_index = Int(floor(r / dr)) + 1
                gr_histogram[bin_index] += 1
            end
        end
    end
    return nothing
end

"""
    normalize_gr!(gr_histogram, num_part, L, n_samples) -> Vector{Float64}

Normalizes an accumulated g(r) histogram in place (ideal-gas normalization
by bin area and sample count) and returns the bin-center radii `r_vals`.
The only value that can't be "returned in place" — everything else is
mutated directly into `gr_histogram`.

Note: this replaces an earlier Stage1-only `normalize_gr` that returned a
fresh (r_vals, gr) tuple instead of mutating. Any call site still using
that non-mutating form needs updating to `r_vals = normalize_gr!(gr, ...)`.
"""
function normalize_gr!(gr_histogram::Vector{Float64}, num_part::Int,
                        L::Float64, n_samples::Int)::Vector{Float64}
    num_bins = length(gr_histogram)
    dr = (L/2) / num_bins
    n = num_part / L^2

    r_vals = [(i - 0.5) * dr for i in 1:num_bins]
    @inbounds for i in 1:num_bins
        gr_histogram[i] /= (n_samples * num_part * n * 2π * r_vals[i] * dr)
    end
    return r_vals
end

"""
    accumulate_density!(n_xy, x_coord, y_coord, L)

In-place accumulation of the real-space 2D density n(x,y) for a single
MC configuration. `n_xy` should be preallocated once outside the MC loop
(zeros(Float64, num_bins, num_bins)) and passed in; call once per sampled
(decorrelated) step.
"""
function accumulate_density!(n_xy::Matrix{Float64},
                              x_coord::Vector{Float64}, y_coord::Vector{Float64}, L::Float64)
    @assert length(x_coord) == length(y_coord) "x_coord and y_coord must have the same length"
    num_bins = size(n_xy, 1)

    @inbounds for i in eachindex(x_coord)
        bin_x = clamp(Int(floor(x_coord[i] / L * num_bins)) + 1, 1, num_bins)
        bin_y = clamp(Int(floor(y_coord[i] / L * num_bins)) + 1, 1, num_bins)
        n_xy[bin_x, bin_y] += 1
    end
    return nothing
end

"""
    density_bin_centers(num_bins, L) -> Vector{Float64}

Bin-center coordinates for the density grid, consistent with the
[0, L) convention used by `accumulate_density!`.
"""
function density_bin_centers(num_bins::Int, L::Float64)::Vector{Float64}
    dx = L / num_bins
    return [(b - 0.5) * dx for b in 1:num_bins]
end

"""
    normalize_density!(n_xy, n_samples, L) -> Vector{Float64}

Normalizes an accumulated density grid in place and returns its bin
centers (see `density_bin_centers`).
"""
function normalize_density!(n_xy::Matrix{Float64}, n_samples::Int, L::Float64)::Vector{Float64}
    num_bins = size(n_xy, 1)
    @assert size(n_xy, 1) == size(n_xy, 2) "n_xy must be square"
    bin_area = (L / num_bins)^2
    norm_factor = 1.0 / (n_samples * bin_area)
    @inbounds for i in eachindex(n_xy)
        n_xy[i] *= norm_factor
    end
    return density_bin_centers(num_bins, L)
end

# TODO:
# angular momentum / current-density estimators for the vortex work.
# Not implemented — dΩ_dx, dΩ_dy, and the trial-wavefunction phase
# gradients (dlnψ_0_dx, dlnψ_0_dy) these depend on don't exist yet.
# Confirm before deleting: still planned, or superseded by something else?
#
# function angular_momentum_z!(Lz_real, Lz_imag, dΩ_dx, dΩ_dy, dlnψ_0_dx, dlnψ_0_dy, x_coord, x0, y_coord, y0)
#     Lz_real = sum((x_coord .- x0).*dΩ_dy - (y_coord .- y0).*dΩ_dx)
#     Lz_imag = sum(-((x_coord .- x0).*dlnψ_0_dy - (y_coord .- y0).*dlnψ_0_dx))
#     return Lz_real, Lz_imag
# end