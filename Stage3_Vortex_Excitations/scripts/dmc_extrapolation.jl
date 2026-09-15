function extrapolate_mixed(O_VMC, O_DMC; guard=1e-12)
    size(O_VMC) == size(O_DMC) || error("The VMC and DMC observables must have the same size.")
    extr_linear = 2 .* O_DMC .- O_VMC
    O_vmc_safe = sign.(O_VMC) .* max.(abs.(O_VMC), guard)
    extr_quadratic = O_DMC .^ 2 ./ O_vmc_safe
    return extr_linear, extr_quadratic
end
    