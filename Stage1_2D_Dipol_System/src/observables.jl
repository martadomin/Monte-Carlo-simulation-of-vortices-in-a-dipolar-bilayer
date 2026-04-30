using Statistics

function blocking_statistics(data::Vector{Float64}, block_size::Int)
    num_blocks = div(length(data), block_size)
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
