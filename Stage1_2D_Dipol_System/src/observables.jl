using Statistics

function blocking_statistics(data::Vector{Float64}, block_size::Int)
    num_blocks = div(length(data), block_size)
    block_means = [mean(data[(i-1)*block_size+1:i*block_size]) for i in 1:num_blocks]
    average = mean(block_means)
    error_average = std(block_means) / sqrt(num_blocks)
    return average, error_average
end
