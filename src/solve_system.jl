
# # block matrix inverse: http://www.cs.nthu.edu.tw/~jang/book/addenda/matinv/matinv/
# function solvesysblock(
#     obj::LeastSquares{T},
#     A_container::GenericA,
#     σ::T,
#     ρ_reciprocal::T,
#     x::Vector{T},
#     y::Vector{T},
#     z::Vector{T},
#     ) where T

#     N = length(x)
#     @assert length(y) == length(z)

#     b1 = σ .* x .- q
#     b2 = z .- ρ_reciprocal .* y

#     A = obj.P + σ .* LinearAlgebra.I
#     B = A_container.matrix'
#     C = A_container.matrix
#     D = diagm(ones(T, N) .* -ρ_reciprocal)
#     inv_D = diagm(ones(T, N) .* ρ_reciprocal)

#     X11 = inv(A-B*inv(D)*C)
#     X12 = -inv(A)*B*inv(D-C*inv(A)*B)
#     X21 = -inv(D)*C*inv(A-B*inv(D)*C)
#     X22 = inv(D-C*inv(A)*B)

#     # out1 = X11*b1 + X12*b2
#     # out2 = X21*b1 + X22*b2

#     # X11*b1
#     term11 = inv(A-B*inv(D)*C)*b1

#     # X12*b2
#     term12 = -inv(A)*B*inv(D-C*inv(A)*B)*b2

#     # X21*b1
#     term21 = -inv(D)*C*inv(A-B*inv(D)*C)*b1

#     # X22*b2
#     term22 = inv(D-C*inv(A)*B)*b2

#     return term11 + term12, term21 + term22
# end


# # block matrix inverse: http://www.cs.nthu.edu.tw/~jang/book/addenda/matinv/matinv/
# function solvesysblock(
#     obj::LeastSquares{T},
#     ::IdentityA,
#     σ::T,
#     ρ_reciprocal,
#     x::Vector{T},
#     y::Vector{T},
#     z::Vector{T},
#     ) where T

#     N = length(x)
#     @assert length(y) == length(z) == N

#     q = obj.q

#     b1 = σ .* x .- q
#     b2 = z .- ρ_reciprocal .* y

#     #A = obj.P + σ .* LinearAlgebra.I
#     #A = obj.PσI
#     inv_A = obj.inv_PσI
#     # B = C = I of size N x N.
#     # D = -ρ_reciprocal .* I
#     inv_D = -ρ .* I

#     # inv(A-inv_D)*b1
#     #tmp1 = (A-inv_D)\b1
#     tmp1 = obj.PσI_minus_inv_D\b1

#     # inv(D-inv(A))*b2
#     #tmp2 = (D-inv_A)\b2
#     tmp2 = obj.D_minus_inv_PσI\b2
    
#     # X11*b1
#     term11 = tmp1

#     # X12*b2
#     term12 = -inv_A*tmp2

#     # X21*b1
#     term21 = -inv_D*tmp1

#     # X22*b2
#     term22 = tmp2

#     return term11 + term12, term21 + term22
# end

# generic version that uses direct solve.
function solvesys(
    obj::LeastSquares{T},
    buffer::DenseSystemBuffer{T},
    σ::T,
    ρ_reciprocal,
    x::Vector{T},
    y::Vector{T},
    z::Vector{T},
    ) where T

    N = length(x)
    @assert length(y) == length(z) == N

    q = obj.q

    b1 = σ .* x .- q
    b2 = z .- ρ_reciprocal .* y
    
    buffer.b[begin:begin+N-1] = b1
    buffer.b[begin+N:end] = b2

    a = buffer.KKT_matrix\buffer.b

    return a[begin:begin+N-1], a[begin+N:end]
end