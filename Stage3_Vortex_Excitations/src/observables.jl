using Statistics
function detect_plateau(block_sizes, sigmas; window_size=4, rtol=0.05)
    """
    Detects the plateau in a blocking analysis.

    A candidate plateau (a flat window of `window_size` points) is only
    accepted if sigma stays within `rtol` of that plateau value for every
    remaining block size, not just within the first flat-looking window.
    This guards against a "shoulder" in sigma(B) — a temporary flat stretch
    from a fast decorrelation mode — being mistaken for true convergence
    before a slower mode pushes sigma up again at larger B.
    """
    if length(sigmas) < window_size
        return block_sizes[end]
    end

    for i in 1:(length(sigmas) - window_size + 1)
        window = sigmas[i:i+window_size-1]
        mean_val = mean(window)
        max_dev = maximum(abs.(window .- mean_val))

        if max_dev / mean_val < rtol
            # Candidate plateau found — verify it holds for ALL later block sizes
            tail = sigmas[i:end]
            tail_mean = mean(tail)
            tail_dev = maximum(abs.(tail .- tail_mean))

            if tail_dev / tail_mean < rtol
                return block_sizes[i]
            end
            # else: false plateau (a shoulder) — keep scanning forward
        end
    end

    println("  [Warning] No clear plateau detected. Using largest block size.")
    return block_sizes[end]
end

function blocking_statistics(data::Vector{Float64}, block_size::Int)
    num_blocks = div(length(data), block_size)
    num_blocks < 2 && return NaN, NaN # Not enough blocks for error estimation
    
    block_means = [mean(data[(i-1)*block_size+1:i*block_size]) for i in 1:num_blocks]
    average = mean(block_means)
    error_average = std(block_means) / sqrt(num_blocks)
    return average, error_average
end

function accumulate_gr!(gr_histogram::Vector{Float64}, x_coord::Vector{Float64},
                        y_coord::Vector{Float64}, L::Float64)
    num_part = length(x_coord)
    num_bins = length(gr_histogram)
    dr = (L/2) / num_bins

    for i in 1:num_part
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

function accumulate_g_AB_r!(gr_histogram::Vector{Float64}, x_coord_A::Vector{Float64},
                          y_coord_A::Vector{Float64}, x_coord_B::Vector{Float64},
                          y_coord_B::Vector{Float64}, L::Float64)
    num_part_A = length(x_coord_A)
    num_part_B = length(x_coord_B)
    num_bins = length(gr_histogram)
    dr = (L/2) / num_bins

    for i in 1:num_part_A
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

function normalize_gr!(gr_histogram::Vector{Float64}, num_part::Int,
                        L::Float64, n_samples::Int)::Vector{Float64}
    num_bins = length(gr_histogram)
    dr = (L/2) / num_bins
    n = num_part / L^2

    r_vals = [(i - 0.5) * dr for i in 1:num_bins]
    @inbounds for i in 1:num_bins
        gr_histogram[i] /= (n_samples * num_part * n * 2π * r_vals[i] * dr)
    end
    return r_vals   # only thing that can't be "returned in place" — still needed by the caller
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
    density_bin_centers(num_bins::Int, L::Float64) -> Vector{Float64}

Bin-center coordinates for the density grid, consistent with the
[0, L) convention used by `accumulate_density!`.
"""
function density_bin_centers(num_bins::Int, L::Float64)::Vector{Float64}
    dx = L / num_bins
    return [(b - 0.5) * dx for b in 1:num_bins]
end

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

# function dΩ_dx()



# end

# function angular_momentum_z!(Lz_real, Lz_imag, dΩ_dx, dΩ_dy, dlnψ_0_dx, dlnψ_0_dy, x_coord, x0, y_coord, y0)
#     Lz_real = sum((x_coord .- x0).*dΩ_dy - (y_coord .- y0).*dΩ_dx)
#     Lz_imag = sum(-((x_coord .- x0).*dlnψ_0_dy - (y_coord .- y0).*dlnψ_0_dx))
#     return Lz_real, Lz_imag
# end

# function current()

# end