
function runOSQPreference(
    P::Matrix{T},
    q::Vector{T},
    lbs::Vector{T},
    ubs::Vector{T},
    A::Matrix{T};
    α = one(T),
    verbose = false,
    ) where T
    #

    # # fresh.
    # Crate OSQP object
    prob2 = OSQP.Model()

    # Setup workspace and change alpha parameter
    P2 = sparse(P)
    G = sparse(A) #sparse(LinearAlgebra.I, N, N)
    OSQP.setup!(prob2; P=P2, q=q, A=G, l=lbs, u=ubs, alpha=α, verbose = verbose)

    # Solve problem
    return OSQP.solve!(prob2)
end