
include("./helpers/reference.jl")

Random.seed!(25)

T = Float64

dic = BSON.load("BLS_prob.bson")
B = dic[:B]
lbs = dic[:lbs]
ubs = dic[:ubs]
y = dic[:observations]



B = randn(4000,4)
y = randn(4000)





N = length(lbs)

P = B'*B
q = -B'*y
A = diagm(ones(T, N))

fill!(ubs, 20.0)

σ = 1e-6
α = 1.6
#σ = 1.0
#α = 1.0
ρ = 1.0
max_iters = 100
residual_tol = 1e-6
x_tol = 1e-6
A_container = OSQPSolver.IdentityA()

config = OSQPSolver.ConfigType(ρ, σ, α, max_iters, residual_tol, x_tol)
buffer = OSQPSolver.DenseSystemBuffer(P, A, σ, ρ)
bounds = OSQPSolver.Box(lbs, ubs)
obj = OSQPSolver.LeastSquares(P, q)
prob = OSQPSolver.Problem(obj, A_container, bounds)

#x0 = (lbs+ubs) ./ 2
x0 = lbs
z0 = A*x0
y0 = zeros(T, N)

sol = OSQPSolver.solve!(
    buffer,
    prob,
    x0,
    z0,
    y0,
    config,
)
@show sol.status, sol.iters_ran
x_star = sol.primal

results_ref = runOSQPreference(P, q, lbs, ubs, A; α = α, verbose = true)
x_OSQP = results_ref.x

f = xx->(dot(xx, P*xx)/2 + dot(q,xx))

@show f(x_star)
@show f(x_OSQP)

####

println("Timing:")
@btime OSQPSolver.solve!(
    buffer,
    prob,
    x0,
    z0,
    y0,
    config,
);
@btime runOSQPreference(P, q, lbs, ubs, A; α = α);


# next, write test script for BLS, and projection problem.
# then work on improving the speed of solve!().