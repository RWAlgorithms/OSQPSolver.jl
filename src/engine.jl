# section 3.2


function solve!(
    prob::Problem,
    x0::Vector{T},
    z0::Vector{T},
    y0::Vector{T},
    config::ConfigType{T},
    ) where T

    max_iters = config.max_iters
    α2 = one(T)-α

    for n = 1:max_iters
        x_tilde_next, v_next = solvesys()

        z_tilde_next = z + 1/ρ .* (v_next-y)
        x_next = α .* x_tilde_next + α2 .* x
        
        tmp = α .* z_tilde_next + α2 .* z + 1/ρ .* y
        z_next = PROJ(tmp, lbs, ubs)

        tmp2 = α .* z_tilde_next + α2 .* z - z_next_next
        y_next = y + ρ .* tmp2
    end

    return nothing
end