
abstract type ConstraintType end

# when A is the identity matrix.
struct BoundConstraints <: ConstraintType end

struct GenericConstraints{MT} <: ConstraintType
    A::MT
end

struct Problem{T <: AbstractFloat, PT <: AbstractMatrix, CT <: ConstraintType}
    y::Vector{T}
    lbs::Vector{T}
    ubs::Vector{T}
    P::PT
    A::CT
end

struct ConfigType{T}
    ρ::T # > 0
    σ::T # > 0 
    α::T # in (0,1)
    max_iters::Int
end

struct Solution{T}
    iters_ran::Int
    primal::Vector{T}
    dual::Vector{T}
    status::Symbol
    minimum::T

    # primal_initial::Vector{T}
    # dual_initial::Vector{T}
end