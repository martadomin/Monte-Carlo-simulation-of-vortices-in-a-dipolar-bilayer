using LinearAlgebra

"""
Detects the plateau in a blocking analysis.

A candidate plateau (a flat window of `window_size` points) is only
accepted if sigma stays within `rtol` of that plateau value for every
remaining block size, not just within the first flat-looking window.
This guards against a "shoulder" in sigma(B) — a temporary flat stretch
from a fast decorrelation mode — being mistaken for true convergence
before a slower mode pushes sigma up again at larger B.
"""
function detect_plateau(block_sizes, sigmas; window_size=4, rtol=0.05)
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

"""
    weighted_extrapolation(x, y, y_err; powers=[1]) -> (intercept, intercept_err, coeffs, fitted_fn)

Weighted least-squares fit of y = a + Σᵢ coeffs[i]·x^(powers[i]), weighted
by 1/y_err², where `powers` fixes exactly which terms are present — the
shape is chosen by the caller from known physics, not inferred from the
data (see examples below).

Returns:
- intercept, intercept_err : the fitted value at x=0 and its uncertainty
                              — the actual extrapolated physical answer
- coeffs                    : the fitted [b, c, ...] coefficients for the
                               non-intercept terms, in the same order as `powers`
- fitted_fn                 : a function of x that evaluates the fitted
                               curve — pass a fine x range to it to draw
                               the fit on a plot

Examples:
- powers=[1]   → y = a + b·x            (linear DMC's O(Δτ) bias;
                                           or a simple 1/num_walkers fit)
- powers=[2]   → y = a + c·x²            (quadratic DMC's O(Δτ²) bias —
                                           NO linear term)
- powers=[1,2] → y = a + b·x + c·x²      (walker-count fit allowing
                                           curvature as a correction)
"""
function weighted_extrapolation(x::Vector{Float64}, y::Vector{Float64},
                                 y_err::Vector{Float64}; powers::Vector{Int}=[1])
    w = 1.0 ./ y_err.^2
    W = Diagonal(w)
    X = hcat(ones(length(x)), (x.^p for p in powers)...)

    β = (X' * W * X) \ (X' * W * y)
    cov = inv(X' * W * X)

    intercept, coeffs = β[1], β[2:end]
    fitted_fn(xv) = intercept + sum(coeffs[i] * xv^powers[i] for i in eachindex(powers)) #defines a function

    return intercept, sqrt(cov[1,1]), coeffs, fitted_fn
end

"""
    choose_largest_consistent(x, y, y_err, intercept, intercept_err; n_sigma=1.0) -> x_chosen

Among tested points (x, y, y_err), returns the largest x whose y is still
consistent with `intercept` within n_sigma combined standard errors —
i.e. |y - intercept| <= n_sigma * sqrt(y_err² + intercept_err²).
Used to pick the cheapest (largest) Δτ still statistically indistinguishable
from the Δτ→0 extrapolated value, rather than picking by eye alone.
"""
function choose_largest_consistent(x::Vector{Float64}, y::Vector{Float64}, y_err::Vector{Float64},
                                    intercept::Float64, intercept_err::Float64; n_sigma::Float64=1.0)
    order = sortperm(x; rev=true) #Gives indexes of x in descending order, so we check largest Δτ first (this is what rev =true does)
    for i in order
        combined_err = sqrt(y_err[i]^2 + intercept_err^2)
        if abs(y[i] - intercept) <= n_sigma * combined_err
            return x[i]
        end
    end
    @warn "No Δτ found consistent with the extrapolated value within $(n_sigma)σ — falling back to the smallest tested Δτ"
    return x[argmin(x)]
end

"""
    combine_runs(values, errors) -> (mean, error)

Weighted mean of several independent measurements of the same quantity,
each with its own statistical error — the correct way to combine repeated
runs (e.g. several save_run entries under one path), not just picking one.
Weight ∝ 1/σ², same inverse-variance weighting as weighted_extrapolation.
"""
function combine_runs(values::Vector{Float64}, errors::Vector{Float64})
    length(values) == 1 && return values[1], errors[1]
    w = 1.0 ./ errors.^2
    mean_val = sum(w .* values) / sum(w)
    combined_err = sqrt(1.0 / sum(w))
    return mean_val, combined_err
end