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
end

function normalize_gr(gr_histogram::Vector{Float64}, num_part::Int,
                      L::Float64, n_samples::Int)::Tuple{Vector{Float64}, Vector{Float64}}
    num_bins = length(gr_histogram)
    dr = (L/2) / num_bins
    n = num_part / L^2

    r_vals = [(i - 0.5) * dr for i in 1:num_bins]
    gr = [gr_histogram[i] / (n_samples * num_part * n * 2π * r_vals[i] * dr)
          for i in 1:num_bins]
    return r_vals, gr
end
