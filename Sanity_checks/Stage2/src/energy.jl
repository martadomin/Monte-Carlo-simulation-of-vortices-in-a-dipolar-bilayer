function energy_estimators(xcoord::Vector{Float64}, ycoord::Vector{Float64},
                            L::Float64, R_match::Float64,
                            Constants::Tuple{Float64, Float64, Float64})::Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64}

    num_part        = length(xcoord)
    drift_x         = zeros(Float64, num_part)
    drift_y         = zeros(Float64, num_part)
    E_kin           = 0.0
    F_drift_sq      = 0.0
    Scalar_term_sum = 0.0
    E_int           = 0.0

    @inbounds for k in 1:num_part
        F_x         = 0.0
        F_y         = 0.0
        scalar_term = 0.0

        for i in 1:num_part
            if i != k
                dx = get_periodic_difference(xcoord[k], xcoord[i], L)
                dy = get_periodic_difference(ycoord[k], ycoord[i], L)
                r  = sqrt(dx^2 + dy^2)

                if r > 1e-10
                    # Compute ONCE per pair
                    du_dr   = u2_first_derivative(r, R_match, L, Constants)
                    d2u_dr2 = u2_second_derivative(r, R_match, L, Constants)

                    F_x         += du_dr * (dx / r)
                    F_y         += du_dr * (dy / r)
                    scalar_term += d2u_dr2 + (du_dr / r)

                    # Interaction only for i < k
                    if i < k && r <= L/2
                        E_int += (dx^2 + dy^2)^(-3/2)
                    end
                end
            end
        end

        drift_x[k]      = F_x
        drift_y[k]      = F_y
        F_drift_sq      += F_x^2 + F_y^2
        Scalar_term_sum += scalar_term
        E_kin           += F_x^2 + F_y^2 + scalar_term
    end

    E_kin_std       = -0.5  * E_kin
    E_kin_drift     =  0.5  * F_drift_sq
    E_kin_laplacian = -0.25 * Scalar_term_sum
    E_total         = E_kin_std + E_int

    return drift_x, drift_y,
           E_total,
           E_kin_drift + E_int,
           E_kin_laplacian + E_int,
           E_kin_std,
           E_int
end