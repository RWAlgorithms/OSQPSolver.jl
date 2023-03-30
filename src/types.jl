
abstract type MatrixAType end

# when A is the identity matrix of size N_vars x N_vars.
struct IdentityA <: MatrixAType end

struct GenericA{MT<:AbstractMatrix} <: MatrixAType
    matrix::MT
end

abstract type BoundsType end

struct UnConstrained <: BoundsType end

struct Box{T} <: BoundsType
    lbs::Vector{T}
    ubs::Vector{T}
end


abstract type ObjectiveType end

# when P is the identity matrix.
struct UnitQuadratic{T} <: ObjectiveType
    q::Vector{T}
end

struct LeastSquares{T} <: ObjectiveType
    #B::MT
    P::Matrix{T}
    q::Vector{T}
end


###############

struct Problem{PT <: ObjectiveType, AT <: MatrixAType, BT <: BoundsType}
    objective::PT
    A::AT
    bounds::BT
end

#####################




####### 

# no sparse matrices.
struct DenseSystemBuffer{T}
    KKT_matrix::Matrix{T}
    b::Vector{T}
end

# see Eq. 24 of (Stellato, 2020)
function DenseSystemBuffer(P::AbstractMatrix{T}, A::AbstractMatrix{T}, σ::T, ρ::T)::DenseSystemBuffer{T} where T
    N = size(P,1)
    @assert size(P,2) == N == size(A,2)

    M = size(A,1)

    buffer = DenseSystemBuffer(Matrix{T}(undef, N+M, N+M), Vector{T}(undef, 2*N))
    resetbuffer!(buffer, P, A, σ, ρ)
    
    return buffer
end

function resetbuffer!(buffer::DenseSystemBuffer{T}, P::AbstractMatrix{T}, A::AbstractMatrix{T}, σ::T, ρ::T) where T

    N = size(P,1)
    @assert size(P,2) == N == size(A,2)
    @assert length(buffer.b) == 2*N

    M = size(A,1)
    @assert size(buffer.KKT_matrix) == (N+M, N+M)

    X11 = P + σ .* LinearAlgebra.I
    X12 = A'
    X21 = A
    X22 = -one(T)/ρ .* diagm(ones(T,N))

    buffer.KKT_matrix[1:N, 1:N] = X11
    buffer.KKT_matrix[1:N, N+1:end] = X12
    buffer.KKT_matrix[N+1:end, 1:N] = X21
    buffer.KKT_matrix[N+1:end, N+1:end] = X22

    return nothing
end

struct BlockSystemBuffer{T}
    # intermediate constants
    PσI::Matrix{T}
    inv_PσI::Matrix{T}

    PσI_minus_inv_D::Matrix{T} # P + σI + ρI, This is A-inv(D) in terms of block matrix inversion.
    D_minus_inv_PσI::Matrix{T} # -1/ρ I - (P + σI). This is D - inv(A) in terms of block matrix inversion.

    invD::Matrix{T} # -ρI.
end

function BlockSystemBuffer(P::Matrix{T}, σ::T, ρ::T)::BlockSystemBuffer{T} where T
    
    inv_D = -ρ .* diagm(ones(T, size(P,1)))
    return BlockSystemBuffer(
        P + σ .* LinearAlgebra.I,
        inv(P + σ .* LinearAlgebra.I),
        P + (σ+ρ) .* LinearAlgebra.I,
        -(P + (σ-one(T)/ρ) .* LinearAlgebra.I),
        inv_D,
    )
end

################

# struct SolverBuffer{T}

# end

#################
struct ConfigType{T}
    ρ::T # > 0
    σ::T # > 0 
    α::T # in (0,1)
    max_iters::Int

    residual_tol::T # > 0.
    x_tol::T
end

struct Solution{T}
    iters_ran::Int
    primal::Vector{T}
    dual::Vector{T}
    status::Symbol
    minimum::T

    primal_residual::Vector{T}
    dual_residual::Vector{T}
    # primal_initial::Vector{T}
    # dual_initial::Vector{T}
end