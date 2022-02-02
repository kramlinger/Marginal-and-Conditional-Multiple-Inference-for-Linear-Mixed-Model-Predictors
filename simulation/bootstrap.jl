
####### SIMULATION OF MARGINAL ###########################################

### +++ INITIALIZATION +++

### import packages
using Random, LinearAlgebra, RDatasets, MixedModels, Distributions, Optim, CSV, Dates, Statistics, CategoricalArrays

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

Q = function(data, y, ybar, x, xbar, n, N, μ0, quantile)

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

    ###     CONDITIONAL

    ### assemble diagonal matrix parts
    l1 = γ^2 / n * σ^2

    l3_1 = 2 * σ^2 * τ^2 / N/ a^3
    l3_2 = 8 / N/ (n - 1) / a^4 / n^2 * τ^2 * σ^2 * (
        + (a^2 - σ^2) / n * (
            + σ^6 / a^2 * n
            - (σ^4  / a^2) * (τ + a / 2 * n)
            - (τ^2 - a / 2 + σ^2 * (1 + n^2)) * (n - 1))
        + σ^10 / a^2 * (n^3 - (n - 1)^2 + 1 / n)
        + a * σ^4 * (n - 1 + 1 / n / 2)  +
        - σ^6 * (2 * n - 1 + 1 / n) +
        + σ^8 / a * (1 + (n-1) * (n^2 + 1) / n) / 2)
    l3_3 = - 8 / N/ n * τ^2 * σ^4 / a^4 * (n - 1) * (a^2 - σ^4)
    l3_4 = - 8 / N/ n * τ^2 * σ^4 / a^3 / (n - 1) * (
        + a * ((n - 1) * n - 1)
        + σ^2 * (n + 1))
    l3 = 2 * (l3_1 + l3_2 + l3_3 + l3_4)

    l4 = 2 * σ^4 / N / a^2 * ( 1 / n + ((n - 1) * τ^2 +
            (σ^2 + τ^2)^2) / ((n - 1) * a^2))

    l5 = 2 * n / (n - 1) / N * σ^4 * (σ^2 - 2 * n * τ^2) / a^2

    ### assemble inverse of covariance matrix with Woodbury
    cΣdiag = maximum([l1 + l3 + l4 - l5, 0.001])
    cΣi = (I - d * inv(cΣdiag * σ^2 *
                inv(vcov1 * (2 * γ * I + n^2 * dd * vcov1)) + dd) * d' ) / cΣdiag


    ###     COMPUTE MEAN
    β = beta(vcov1, x, xbar, y, ybar, σ, τ, n, γ) # = fixef(model)
    μ = (1 - γ) * (β[1] .+ β[2] .* xbar) + γ * ybar

    ###     SAVE PIVOTS
    MarginalPivot = sum((μ - μ0)' .* ((μ - μ0)' * Σi), dims = 2)[1]
    ConditionalPivot = sum((μ - μ0)' .* ((μ - μ0)' * cΣi), dims = 2)[1]

    if quantile == true

        ###     NONCENTRALITY PARAMETER
        ay = d * β - (1 - γ) * ybar
        λ = maximum([(ay' * cΣi * ay)[1] - tr(cΣi) * σ^2 / a^2 / n, 0])

        ###     QUANTILES
        qc = cquantile(NoncentralChisq(m, λ), 0.05) # quantile of alpha = 0.05

        return [MarginalPivot ConditionalPivot qc β[1] β[2] τ σ Σdiag cΣdiag]
    else
        return [MarginalPivot ConditionalPivot]
    end
end


### fix parameters
p = 2
nsim = 5000
R = 300
Random.seed!(2021)

### +++ PER SPECIFICATION SETTING +++

### to compile
τ0 = 8^0.5 
σ0 = 4^0.5
m = 50
n = 10
REdistr = "normal"
Edistr = "normal"

### variable parameters
N = m * n

β0 = [-4.9, 0.03]
x = 10 * randn(N)
xbar = vec(sum(reshape(x, n, :), dims=1) / n)

qm = cquantile(Chisq(m), 0.05) # quantile of alpha = 0.05

### conditional evaluation
if REdistr == "normal"
    v = τ0 * randn(m)
elseif REdistr == "tdist"
    v = rand(TDist(2), m)
elseif REdistr == "chi2"
    v = rand(Chisq(τ0^2 / 2), m) .- τ0^2 / 2
end
u = β0[1] .+ β0[2] * x + repeat(v, inner = n)
μ0 = β0[1] .+ β0[2] .* xbar + v

### obtain conditional variance
σ1 = zeros(10000, 1)
τ1 = zeros(10000, 1)
for i in 1:10000
    if Edistr == "normal"
        e = σ0 * randn(N)
    elseif Edistr == "tdist"
        e = rand(TDist(2), N)
    elseif Edistr == "chi2"
        e = rand(Chisq(σ0^2 / 2), N) .- σ0^2 / 2
    end
    y = u + e # conditional
    data = DataFrame(x = x, v = repeat(categorical(v), inner = n), y = y)
    model = fit(MixedModel, @formula(y ~ x + (1 | v)), data)
    τ1[i] = VarCorr(model).σρ.v.σ[1]
    σ1[i] = VarCorr(model).s
end
σ1 = mean(σ1) # should be (close to) σ0, but theory for conditional convergence is missing
τ1 = mean(τ1) # should be (close to) τ0, but theory for conditional convergence is missing

### simulation
Coverage = zeros(6)
λ = zeros(nsim)
Volume = ones(nsim, 6)
QM1 = zeros(nsim)
QC1 = zeros(nsim)
QM2 = zeros(nsim)
QC2 = zeros(nsim)
for i in 1:nsim

    ### +++ EMPIRICAL +++

    ### variable data
    if Edistr == "normal"
        e = σ0 * randn(N)
    elseif Edistr == "tdist"
        e = rand(TDist(2), N)
    elseif Edistr == "chi2"
        e = rand(Chisq(σ0^2 / 2), N) .- σ0^2 / 2
    end
    y = u + e # conditional
    ybar = vec(sum(reshape(y, n, :), dims=1) / n)
    data = DataFrame(x = x, v = repeat(categorical(v), inner = n), y = y)

    MP_bt_Alg1 = zeros(R-1)
    CP_bt_Alg1 = zeros(R-1)

    β1 = 1;  β2 = 1; τ = 1;  σ = 1

    MP, CP, qc, β1, β2, τ, σ, Σdiag, cΣdiag= Q(data, y, ybar, x, xbar, n, N, μ0, true)

    for r in 2:R

        v_bt = τ * randn(m)
        u_bt = β1 .+ β2 * x + repeat(v_bt, inner = n)
        μ0_bt = β1 .+ β2 .* xbar + v_bt

            ### generate new bootstrap data
            y_bt = u_bt + σ * randn(N) # conditional
            ybar_bt = vec(sum(reshape(y_bt, n, :), dims = 1) / n)
            data_bt = DataFrame(x = x, v = repeat(categorical(v_bt), inner = n), y = y_bt)

            MP_bt_Alg1[r-1], CP_bt_Alg1[r-1] = Q(data_bt, y_bt, ybar_bt, x, xbar, n, N, μ0_bt, false)
    end


    MP_bt_Alg2 = zeros(R-1)
    CP_bt_Alg2 = zeros(R-1)

    v_bt = τ * randn(m)
    u_bt = β1 .+ β2 * x + repeat(v_bt, inner = n)
    μ0_bt = β1 .+ β2 .* xbar + v_bt

    for r in 2:R

            ### generate new bootstrap data
            y_bt = u_bt + σ * randn(N) # conditional
            ybar_bt = vec(sum(reshape(y_bt, n, :), dims = 1) / n)
            data_bt = DataFrame(x = x, v = repeat(categorical(v_bt), inner = n), y = y_bt)

            MP_bt_Alg2[r-1], CP_bt_Alg2[r-1] = Q(data_bt, y_bt, ybar_bt, x, xbar, n, N, μ0_bt, false)
    end

    ###     COMPUTE BOOTSTRAP QUANTILES
    QM1[i] = quantile!(MP_bt_Alg1, 0.95)
    QC1[i] = quantile!(CP_bt_Alg1, 0.95)


    ###     COMPUTE BOOTSTRAP QUANTILES
    QM2[i] = quantile!(MP_bt_Alg2, 0.95)
    QC2[i] = quantile!(CP_bt_Alg2, 0.95)

    ###     COMPUTE COVERAGE
    Coverage[1] += (MP .< qm) / nsim
    Coverage[2] += (CP .< qc) / nsim
    Coverage[3] += (MP .< QM1[i]) / nsim
    Coverage[4] += (CP .< QC1[i]) / nsim
    Coverage[5] += (MP .< QM2[i]) / nsim
    Coverage[6] += (CP .< QC2[i]) / nsim

    ###     COMPUTE RELATIVE VOLUME
    # (via approximations as the determinant of large matrices is numerically unstable)
    Volume[i, 2] = (Σdiag / cΣdiag * qm / qc)^(-m / 2)
    Volume[i, 3] = (qm / QM1[i])^(-m / 2)
    Volume[i, 4] = (Σdiag / cΣdiag * qc / QC1[i])^(-m / 2)
    Volume[i, 5] = (qm / QM2[i])^(-m / 2)
    Volume[i, 6] = (Σdiag / cΣdiag * qc / QC2[i])^(-m / 2)

end
