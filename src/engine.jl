# section 3.2


function solve!(
    buffer::DenseSystemBuffer{T},
    prob::Problem,
    x0::Vector{T},
    z0::Vector{T},
    y0::Vector{T},
    config::ConfigType{T},
    ) where T

    N = length(x0)
    @assert length(z0) == length(y0) == N

    max_iters = config.max_iters
    objective = prob.objective
    bounds = prob.bounds
    A = prob.A
    
    # setup constants.
    α, ρ, σ = config.α, config.ρ, config.σ
    residual_tol, x_tol = config.residual_tol, config.x_tol
    α2 = one(T)-α
    ρ_reciprocal = one(T)/ρ

    # variables.
    x_next = Vector{T}(undef, N)
    y_next = Vector{T}(undef, N)
    z_next = Vector{T}(undef, N)

    x = copy(x0)
    y = copy(y0)
    z = copy(z0)

    # # intermediate buffers.
    # α_z_tilde_next = Vector{T}(undef, N)
    # α2_z = Vector{T}(undef, N)

    for n = 1:max_iters
        # # stopping condition check.
        r_primal, r_dual = getresiduals(A, objective, x, y, z)
        
        if norm(r_primal,Inf) + norm(r_dual,Inf) <  residual_tol

            f_x = evalobjective(objective, x)
            return Solution(n, x, y, :converged, f_x, r_primal, r_dual)
        end

        # # update.
        x_tilde_next, v_next = solvesys(objective, buffer, σ, ρ_reciprocal, x, y, z)

        z_tilde_next = z + ρ_reciprocal .* (v_next-y)
        # α_z_tilde_next[:] = z + ρ_reciprocal .* (v_next .- y)
        # α_z_tilde_next[:] = α .* α_z_tilde_next
        # α2_z[:] = α2 .* z

        x_next = α .* x_tilde_next + α2 .* x
        
        z_next = α .* z_tilde_next + α2 .* z + ρ_reciprocal .* y
        #z_next = α_z_tilde_next + α2_z + ρ_reciprocal .* y
        project!(z_next, bounds)

        tmp2 = α .* z_tilde_next + α2 .* z - z_next
        y_next = y + ρ .* tmp2

        # # book keep.
        x, x_next = x_next, x
        y, y_next = y_next, y
        z, z_next = z_next, z

        if norm(x-x_next,Inf) < x_tol
            f_x = evalobjective(objective, x)
            return Solution(n, x, y, :iterate_tol_reached, f_x, r_primal, r_dual)
        end
    end

    r_primal, r_dual = getresiduals(A, objective, x, y, z)
    f_x = evalobjective(objective, x)
    return Solution(max_iters, x, y, :max_iters_reached, f_x, r_primal, r_dual)
end


function project!(x::Vector{T}, C::Box{T}) where T <: Real
    for d in eachindex(x)
        x[d] = clamp(x[d], C.lbs[d], C.ubs[d])
    end
    return nothing
end

function project!(x, ::UnConstrained)
    return nothing
end

function getresiduals(
    ::IdentityA,
    obj::LeastSquares{T},
    x::Vector{T},
    y::Vector{T},
    z::Vector{T},
    ) where T

    r_primal = x-z
    r_dual = obj.P*x .+ obj.q + y

    return r_primal, r_dual
end

function getresiduals(
    A_container::GenericA{T},
    obj::LeastSquares{T},
    x::Vector{T},
    y::Vector{T},
    z::Vector{T},
    ) where T

    A = A_container.matrix

    r_primal = A*x-z
    r_dual = obj.P*x .+ obj.q + A'*y

    return r_primal, r_dual
end

function evalobjective(obj::LeastSquares{T}, x::Vector{T}) where T
    P = obj.P
    q = obj.q
    return dot(x,P*x)/2 + dot(q,x)
end

function evalobjective(obj::UnitQuadratic{T}, x::Vector{T}) where T
    q = obj.q
    return dot(x,x)/2 + dot(q,x)
end