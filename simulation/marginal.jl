
####### SIMULATION OF MARGINAL ###########################################

### +++ INITIALIZATION +++

### import packages
using Random, LinearAlgebra, RDatasets, MixedModels, Distributions, Optim, CSV, Dates

### define functions
vcov2 = function(x, xbar, σ, τ, n, γ)
    m = length(xbar)
    a = γ / τ^2 * m # 1-1 entry
    d = (sum(x.^2) - n * γ * sum(xbar.^2) ) / σ^2 # 2-2 entry
    b = γ / τ^2 / n * sum(x) # 2-1 entry
    [d -b; -b a] / (a * d - b^2) # return inverse (X^t V^-1 X)^-1
end

beta = function(vcov1, x, xbar, y, ybar, σ, τ, n, γ)
    c1 = γ / τ^2 * sum(ybar) # = 1^t V^-1 y
    c2 = (x' * y - n * γ * sum(ybar .* xbar))/σ^2 # = x^t V^-1 y
    vcov1 * [c1 c2]'  # (X^tV^-1X)^-1 X^tV^-1y
end

### fix parameters
p = 2
nsim = 10000
Δ = [0:0.1:2;]
Random.seed!(12345)

### +++ PER SPECIFICATION SETTING +++

### to compile
τ1 = 8^0.5
σ1 = 4^0.5
m = 10
n = 1

### variable parameters
N = m * n

β0 = [-4.9, 0.03] # rand(Uniform(0, 1), p)
x = rand(Uniform(0, 10), N)
xbar = vec(sum(reshape(x, n, :), dims=1) / n)

### +++ KNOWN DELTA +++

###     MARGINAL COVARIANCE MATRIX

### helpful parameters
γ0 = τ1^2 / (τ1^2 + σ1^2 / n)
vcov0 = vcov2(x, xbar, σ1, τ1, n, γ0) # = (X^tV^-1X)^-1
d0 = (1 - γ0) .* hcat(fill(1, m), xbar) # d' = l' - b'X
dd0 = d0' * d0

### assemble covariance matrix parts
g1 = γ0 * σ1^2 / n

### assemble inverse of covariance matrix with Woodbury
Σ0diag = g1
Σ0i = (I - d0 * inv(Σ0diag * inv(vcov0) + dd0) * d0' ) / Σ0diag

###     QUANTILES
qm = cquantile(Chisq(m), 0.05) # quantile of alpha = 0.05

### simulation
Coverage = zeros(2, length(Δ))
for i in 1:nsim

    ### +++ EMPIRICAL +++

    ### marginal evaluation
    v = τ1 * randn(m)
    u = β0[1] .+ β0[2] * x + repeat(v, inner = n)
    μ0 = β0[1] .+ β0[2] .* xbar + v

    ### variable data
    y = u + σ1 * randn(N) # conditional
    ybar = vec(sum(reshape(y, n, :), dims=1) / n)
    data = DataFrame(x = x, v = repeat(categorical(v), inner = n), y = y)

    ### fit model
    model = fit(MixedModel, @formula(y ~ x + (1 | v)), data)

    ### make sure the variance parameters are positive
    τ = maximum([VarCorr(model).σρ.v.σ[1], 0.001])
    σ = maximum([VarCorr(model).s, 0.001])

    ###     MARGINAL

    ### helpful parameters
    γ = τ^2 / (τ^2 + σ^2 / n)
    vcov1 = vcov2(x, xbar, σ, τ, n, γ) # = (X^tV^-1X)^-1
    d = (1-γ) .* hcat(fill(1, m), xbar) # d' = l' - b'X
    dd = d' * d
    a = n * τ^2 + σ^2

    ### assemble diagonal matrix parts
    g1 = γ * σ^2 / n
    g3 = 2 * σ^4 / τ^2 * γ / N/ (n - 1)

    ### assemble inverse of covariance matrix with Woodbury
    Σdiag = g1 + 2 * g3
    #Σ = Σdiag * I + d * vcov1 * d'
    Σi = (I - d * inv(Σdiag * inv(vcov1) + dd) * d' ) / Σdiag

    ###     COMPUTE MEAN
    β1 = beta(vcov0, x, xbar, y, ybar, σ1, τ1, n, γ0)
    μ1 = (1 - γ0) * (β1[1] .+ β1[2] .* xbar) + γ0 * ybar

    β2 = beta(vcov1, x, xbar, y, ybar, σ, τ, n, γ) # = fixef(model)
    μ2 = (1 - γ) * (β2[1] .+ β2[2] .* xbar) + γ * ybar

    ###     COMPUTE STATISTICS
    t1 = (μ1 - μ0) .+ Δ'
    t2 = (μ2 - μ0) .+ Δ'
    Coverage[1,:] += (sum(t1' .* (t1' * Σ0i), dims = 2) .< qm) / nsim
    Coverage[2,:] += (sum(t2' .* (t2' * Σi), dims = 2) .< qm) / nsim
end
