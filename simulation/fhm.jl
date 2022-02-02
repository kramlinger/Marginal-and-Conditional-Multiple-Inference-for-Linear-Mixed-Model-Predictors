
####### SIMULATION OF MARGINAL ###########################################

### +++ INITIALIZATION +++

### import packages
using Random, LinearAlgebra, RDatasets, MixedModels, Distributions, Optim, CSV, Dates, CategoricalArrays, GLM


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


###     SELECT PARAMETERS
m = 10
n = 1
τ0 = 8^0.5
σ0 = 4^0.5

### fix parameters
p = 2
nsim = 10000
Δ = 0
Random.seed!(123)

### +++ PER SPECIFICATION SETTING +++
REdistr = "normal"

### variable parameters
N = m * n

β0 = [-4.9, 0.03]
x = 10 * randn(N)
xbar = vec(sum(reshape(x, n, :), dims=1) / n)
X = hcat(fill(1, N), x)

### conditional evaluation
if REdistr == "normal"
    v = τ0 * randn(m)
elseif REdistr == "tdist"
    v = rand(TDist(2), m)
elseif REdistr == "chi2"
    v = rand(Chisq(τ0^2 / 2), m) .- τ0^2 / 2 - m
end
u = β0[1] .+ β0[2] * x + repeat(v, inner = n)
μ0 = β0[1] .+ β0[2] .* xbar + v

hatmatrix = inv(X' * X) * X'
detX = det(X' * X)
loglik = function(α, y)
    ξ = exp(α[1])
    beta = hatmatrix * y
    residual = y - X * beta
    (p * log(detX / ξ) + N * log(ξ) + residual' * residual / ξ)
end

τ1 = zeros(10000, 1)
for i in 1:10000
    y = u + σ0 * randn(N)
    τ1[i] = sqrt(maximum([exp(optimize(α -> loglik(α,y), [1.0], BFGS()).minimizer[1]) - σ0^2 0]))
end
τ1 = mean(τ1) # should be (close to) τ0, but theory for conditional convergence is missing
σ1 = σ0  #mean(σ1)

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

###     CONDITIONAL COVARIANCE MATRIX
l1 = γ0^2 / n * σ1^2

### assemble inverse of covariance matrix with Woodbury
cΣ0diag = l1
cΣ0i = (I - d0 * inv(cΣ0diag * σ1^2 *
            inv(vcov0 * (2 * γ0 * I + n^2 * dd0 * vcov0)) + dd0) * d0' ) / cΣ0diag

###     NONCENTRALITY PARAMETER
AZv = (γ0 - 1) .* v + n / τ1^2 .* d0 * vcov0 * d0' * v
λ0 = AZv' * cΣ0i * AZv

###     QUANTILES
qc0 = cquantile(NoncentralChisq(m, λ0), 0.05) # quantile of alpha = 0.05
qm = cquantile(Chisq(m), 0.05) # quantile of alpha = 0.05

### simulation
Coverage = zeros(5, length(Δ))
λ = zeros(nsim, 2)
Volume = ones(nsim, 5)
Τ = zeros(nsim, 1)
Statistic = zeros(nsim, 4)

for i in 1:nsim

    # i = 1
    ### +++ EMPIRICAL +++
    y = u + σ0 * randn(N)
    τ = sqrt(maximum([exp(optimize(α -> loglik(α,y), [1.0], BFGS()).minimizer[1]) - σ0^2 0.01]))
    Τ[i,1] = τ

    ybar = vec(sum(reshape(y, n, :), dims=1) / n)

    ### make sure the variance parameters are positive
    σ = σ1

    ###     MARGINAL

    ### helpful parameters
    γ = τ^2 / (τ^2 + σ^2 / n)
    vcov1 = vcov2(x, xbar, σ, τ, n, γ) # = (X^tV^-1X)^-1
    d = (1-γ) .* hcat(fill(1, m), xbar) # d' = l' - b'X
    dd = d' * d
    a = n * τ^2 + σ^2

    ### assemble diagonal matrix parts
    g1 = γ * σ^2 / n
    g3 = 2 * σ^4 / τ^2 * γ / N

    ### assemble inverse of covariance matrix with Woodbury
    Σdiag = g1 + 2 * g3
    Σ = Σdiag * I + d * vcov1 * d'
    Σi = (I - d * inv(Σdiag * inv(vcov1) + dd) * d' ) / Σdiag

    ###     CONDITIONAL

    ### assemble diagonal matrix parts
    l1 = γ^2 / n * σ^2

    l3_1 = 2 * σ^2 * τ^2 / N/ a^3
    l3_2 = 8 / N / a^4 / n^2 * τ^2 * σ^2 * (
        + (a^2 - σ^2) / n * (
            + σ^6 / a^2 * n
            - (σ^4  / a^2) * (τ + a / 2 * n))
        + σ^10 / a^2 * (n^3 + 1 / n)
        + a * σ^4 * (n - 1 + 1 / n / 2)  +
        - σ^6 * (2 * n - 1 + 1 / n) +
        + σ^8 / a * (1 + (n-1) * (n^2 + 1) / n) / 2)
    l3_3 = 0
    l3_4 = - 8 / N/ n * τ^2 * σ^4 / a^3 * (
        - a
        + σ^2 * (n + 1))
    l3 = 2 * (l3_1 + l3_2 + l3_3 + l3_4)

    l4 = 2 * σ^4 / N / a^2 * ( 1 / n + 1)

    l5 = 2 * n / N * σ^4 * (σ^2 - 2 * n * τ^2) / a^2

    ### assemble inverse of covariance matrix with Woodbury
    cΣdiag = maximum([l1 + l3 + l4 - l5, 1/100000])
    cΣi = (I - d * inv(cΣdiag * σ^2 *
                inv(vcov1 * (2 * γ * I + n^2 * dd * vcov1)) + dd) * d' ) / cΣdiag


    ###     COMPUTE MEAN
    β1 = beta(vcov0, x, xbar, y, ybar, σ1, τ1, n, γ0)
    μ1 = (1 - γ0) * (β1[1] .+ β1[2] .* xbar) + γ0 * ybar

    β2 = beta(vcov1, x, xbar, y, ybar, σ, τ, n, γ) # = fixef(model)
    μ2 = (1 - γ) * (β2[1] .+ β2[2] .* xbar) + γ * ybar

    ###     NONCENTRALITY PARAMETER
    ay = d0 * β1 - (1 - γ0) * ybar
    λ[i, 1] = maximum([(ay' * cΣ0i * ay)[1] - tr(cΣ0i) * γ0^2 / n^3 / τ1^4 * σ1^2, 0])

    ay = d * β2 - (1 - γ) * ybar
    λ[i, 2] = maximum([(ay' * cΣi * ay)[1] - tr(cΣi) * σ^2 / a^2 / n, 0])



    ###     QUANTILES
    qc1 = cquantile(NoncentralChisq(m, λ[i, 1]), 0.05) # quantile of alpha = 0.05
    qc2 = cquantile(NoncentralChisq(m, λ[i, 2]), 0.05) # quantile of alpha = 0.05

    ###     COMPUTE STATISTICS
    t1 = (μ1 - μ0) .+ Δ'
    t2 = (μ2 - μ0) .+ Δ'

    Statistic[i,:] = [sum(t1' .* (t1' * Σ0i), dims = 2)[1],  #marginal known
                      sum(t2' .* (t2' * Σi), dims = 2)[1], #marginal REML
                      sum(t1' .* (t1' * cΣ0i), dims = 2)[1], #conditional known
                      sum(t2' .* (t2' * cΣi), dims = 2)[1]] #conditional REML

    Coverage[1,1] += (Statistic[i,1] .< qm) / nsim
    Coverage[2,1] += (Statistic[i,2] .< qm) / nsim
    Coverage[3,1] += (Statistic[i,3] .< qc0) / nsim
    Coverage[4,1] += (Statistic[i,3] .< qc1) / nsim
    Coverage[5,1] += (Statistic[i,4] .< qc2) / nsim

    ###     COMPUTE RELATIVE VOLUME
    # (via approximations as the determinant of large matrices is numerically unstable)
    Volume[i, 3] = (Σ0diag / cΣ0diag * qm / qc0)^(-m / 2)
    Volume[i, 4] = (Σ0diag / cΣ0diag * qm / qc1)^(-m / 2)
    Volume[i, 5] = (Σdiag / cΣdiag * qm / qc2)^(-m / 2)
end
