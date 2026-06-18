using DelimitedFiles, Plots, LaTeXStrings, PGFPlotsX

function RK4(f::Function, r0::Float64, y0::Vector{Float64}, Δ::Float64, h::Float64, energy_b::Float64)::Vector{Float64}
    k1 = f(r0, y0, h, energy_b)
    k2 = f(r0 + Δ/2, y0 + Δ/2 * k1, h, energy_b)
    k3 = f(r0 + Δ/2, y0 + Δ/2 * k2, h, energy_b)
    k4 = f(r0 + Δ, y0 + Δ * k3, h, energy_b)
    return y0 + (Δ/6) * (k1 + 2*k2 + 2*k3 + k4)
end

function V_AB(r::Float64, h::Float64)::Float64
    r2 = r^2
    h2 = h^2
    return (r2 - 2h2) / (r2 + h2)^(5/2)
end

function ode_rhs(r::Float64, y::Vector{Float64}, h::Float64, energy_b::Float64)::Vector{Float64}
    dy = zeros(2)
    dy[1] = y[2]
    dy[2] = -y[2]/r + (V_AB(r, h) - energy_b) * y[1]
    return dy
end

function initial_conditions(r0::Float64, h::Float64, energy_b::Float64)::Vector{Float64}
    VAB = V_AB(0.0, h)
    y1 = 1.0 + (VAB - energy_b) * r0^2 / 4
    y2 = (VAB - energy_b) * r0 / 2
    return [y1, y2]
end

function shooting_method(r0::Float64, h::Float64, energy_b::Float64, r_max::Float64, Δ::Float64)::Float64
    y0 = initial_conditions(r0, h, energy_b)
    r = r0
    while r < r_max
        y0 = RK4(ode_rhs, r, y0, Δ, h, energy_b)
        r += Δ
    end
    return y0[2]
end

function find_energy_b(r0::Float64, h::Float64, r_max::Float64,
                       Δ::Float64, tol::Float64)::Float64

    F(eb) = shooting_method(r0, h, eb, r_max, Δ)

    E_lo = -2.0 / h^3
    E_hi = -1e-14

    F_lo = F(E_lo)
    F_hi = F(E_hi)

    if F_lo * F_hi > 0
        error("No sign change in bracket (-2/h³, 0) for h=$h, R0=$r_max. F_lo=$F_lo, F_hi=$F_hi")
    end

    while (E_hi - E_lo) > tol
        E_mid = (E_lo + E_hi) / 2.0
        F_mid = F(E_mid)
        if F_lo * F_mid < 0
            E_hi = E_mid
        else
            E_lo = E_mid
            F_lo = F_mid
        end
    end

    return (E_lo + E_hi) / 2.0
end

function build_fAB(h::Float64, R0::Float64, r_min::Float64,
                   Δ::Float64, tol::Float64, nr0sq::Float64, N::Int)

    Δ = min(Δ, h^(3/2) / 20.0)

    energy_b = find_energy_b(r_min, h, R0, Δ, tol)

    n_steps = round(Int, (R0 - r_min) / Δ)
    r_grid  = Vector{Float64}(undef, n_steps + 1)
    psi     = Vector{Float64}(undef, n_steps + 1)
    dpsi    = Vector{Float64}(undef, n_steps + 1)

    y0        = initial_conditions(r_min, h, energy_b)
    r_grid[1] = r_min
    psi[1]    = y0[1]
    dpsi[1]   = y0[2]

    r = r_min
    for i in 1:n_steps
        y0          = RK4(ode_rhs, r, y0, Δ, h, energy_b)
        r          += Δ
        r_grid[i+1] = r
        psi[i+1]    = y0[1]
        dpsi[i+1]   = y0[2]
    end

    # Return nothing instead of erroring — caller decides what to do
    if psi[1] < 0.0 || any(i -> psi[i] * psi[i-1] < 0.0, 2:length(psi))
        return nothing
    end

    psi_R0 = psi[end]
    psi  ./= psi_R0
    dpsi ./= psi_R0

    u_prime       = Vector{Float64}(undef, length(r_grid))
    u_doubleprime = Vector{Float64}(undef, length(r_grid))
    for i in 1:length(r_grid)
        u_prime[i]       = dpsi[i] / psi[i]
        d2psi_i          = (V_AB(r_grid[i], h) - energy_b) * psi[i] - dpsi[i]/r_grid[i]
        u_doubleprime[i] = d2psi_i / psi[i] - u_prime[i]^2
    end

    return r_grid, psi, u_prime, u_doubleprime, energy_b
end