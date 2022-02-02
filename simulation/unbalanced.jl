
####### SIMULATION OF MARGINAL ###########################################

### +++ INITIALIZATION +++

### import packages
using Random, LinearAlgebra, RDatasets, MixedModels, Distributions, Optim, CSV, Dates

### define functions
vcov2 = function(x, xbar, σ, τ, nj, mj, γ)
    m = length(xbar)
    a = sum(γ / τ^2 .* mj)  # 1-1 entry
    d = (sum(x.^2) - sum(rep(γ .* nj, mj) .* xbar.^2)) / σ^2 # 2-2 entry
    b = sum(rep(γ, mj) .* xbar) / τ^2 # 2-1 entry
    [d -b; -b a] / (a * d - b^2) # return inverse (X^t V^-1 X)^-1
end

beta = function(vcov1, x, xbar, y, ybar, σ, τ, nj, mj, γ)
    c1 = sum(rep(γ, mj) .* ybar) / τ^2 # = 1^t V^-1 y
    c2 = (x' * y - sum(rep(nj .* γ, mj) .* ybar .* xbar)) / σ^2 # = x^t V^-1 y
    vcov1 * [c1 c2]'  # (X^tV^-1X)^-1 X^tV^-1y
end

rep = function(x, y)
    n = length(y)
    out = []
    for i in 1:n
        append!(out, fill(x[i], y[i]))
    end
    out
end

max1 = function(x, y)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i] = maximum([x[i], y])
    end
    out
end

τlist = repeat([4, 6, 8].^0.5, inner = 3)
σlist = repeat([8, 6, 4].^0.5, inner = 3)

mjlist = repeat([[10, 30, 10],[20, 20, 10],[30, 10, 10]], 3)

for j in 1:length(mjlist)

    ###     SELECT PARAMETERS
    mj = mjlist[j]
    τ0 = τlist[j]
    σ0 = σlist[j]

    ### fix parameters
    p = 2
    nsim = 10000
    Δ = [0:0.1:2;]
    Random.seed!(999)

    ### +++ PER SPECIFICATION SETTING +++

    ### to compile
    if !@isdefined(τ0) τ0 = 1 end
    if !@isdefined(σ0) σ0 = 1 end
    if !@isdefined(mj) mj = 1 end

    ### variable parameters
    nj = [1, 10, 50]
    m = sum(mj)
    N = sum(nj .* mj)


    β0 = [-4.9, 0.03] # rand(Uniform(0, 1), p)
    x = 10 * randn(N)
    xbar = by(DataFrame(x = x, v = rep(1:m, rep(nj, mj))), :v, :x => mean)[:,2]

    ### conditional evaluation
    v = τ0 * randn(m)
    u = β0[1] .+ β0[2] * x + rep(v, rep(nj, mj))
    μ0 = β0[1] .+ β0[2] .* xbar + v

    ### obtain conditional variance
    σ1 = zeros(10000, 1)
    τ1 = zeros(10000, 1)
    for i in 1:10000
        data = DataFrame(x = x, v = rep(v, rep(nj, mj)), y = u + σ0 * randn(N))
        model = fit(MixedModel, @formula(y ~ x + (1 | v)), data)
        τ1[i] = VarCorr(model).σρ.v.σ[1]
        σ1[i] = VarCorr(model).s
    end
    σ1 = mean(σ1) # should be (close to) σ0, but theory for conditional convergence is missing
    τ1 = mean(τ1) # should be (close to) τ0, but theory for conditional convergence is missing

    ### +++ KNOWN DELTA +++

    ###     MARGINAL COVARIANCE MATRIX

    ### helpful parameters
    γ0 = τ1^2 ./ (τ1^2 .+ σ1^2 ./ nj)
    vcov0 = vcov2(x, xbar, σ1, τ1, nj, mj, γ0) # = (X^tV^-1X)^-1
    d0 = rep(1 .- γ0, mj) .* hcat(fill(1, m), xbar) # d' = l' - b'X
    dd0 = d0' * d0

    ### assemble covariance matrix parts
    g1 = rep(γ0 .* σ1^2 ./ nj, mj)

    ### assemble inverse of covariance matrix with Woodbury
    Σ0diag = convert(Array{Float64,1}, g1)
    Σ0i = (I - (d0 ./ Σ0diag) *
                inv(inv(vcov0) +  d0' * (d0 ./ Σ0diag)) * d0') * Diagonal(Σ0diag.^-1)


    ###     CONDITIONAL COVARIANCE MATRIX
    l1 = rep(γ0.^2 ./ nj * σ1^2, mj)

    ### assemble inverse of covariance matrix with Woodbury
    cΣ0diag = convert(Array{Float64,1}, l1)
    cΣ0 = Diagonal(cΣ0diag) + d0 * (vcov0 + vcov0 * [1 sum(xbar); sum(xbar) sum(xbar.^2)] *
                vcov0 .* sum(nj .* γ0 .* (1 .- γ0) .* mj) ./ σ1^2) * d0' + (rep(γ0, mj) .* d0) *
                vcov0 * d0' + ((rep(γ0, mj) .* d0) * vcov0 * d0')'
    cΣ0i = inv(cΣ0)

    ###     NONCENTRALITY PARAMETER
    AZv = rep(γ0 .- 1, mj) .* v + rep(nj, mj) .* d0 * vcov0 * d0' * v ./ τ1^2
    λ0 = AZv' * cΣ0i * AZv

    ###     QUANTILES
    qc0 = cquantile(NoncentralChisq(m, λ0), 0.05) # quantile of alpha = 0.05
    qm = cquantile(Chisq(m), 0.05) # quantile of alpha = 0.05

    ### simulation
    Coverage = zeros(5, length(Δ))
    λ = zeros(nsim, 2)
    Volume = ones(nsim, 5)
    for i in 1:nsim

        ### +++ EMPIRICAL +++

        ### variable data
        y = u + σ0 * randn(N) # conditional
        ybar = by(DataFrame(y = y, v = rep(1:m, rep(nj, mj))), :v, :y => mean)[:,2]
        data = DataFrame(x = x, v = rep(v, rep(nj, mj)), y = u + σ0 * randn(N))

        ### fit model
        model = fit(MixedModel, @formula(y ~ x + (1 | v)), data)

        ### make sure the variance parameters are positive
        τ = maximum([VarCorr(model).σρ.v.σ[1], 0.001])
        σ = maximum([VarCorr(model).s, 0.001])

        ###     MARGINAL

        ### helpful parameters
        γ = τ^2 ./ (τ^2 .+ σ^2 ./ nj)
        vcov1 = vcov2(x, xbar, σ, τ, nj, mj, γ) # = (X^tV^-1X)^-1
        d = rep(1 .- γ, mj)  .* hcat(fill(1, m), xbar) # d' = l' - b'X
        #dd = d' * d
        aj = nj .* τ^2 .+ σ^2

        Ivv = sum(rep(nj.^2 ./ aj.^2, mj)) / 2
        Ive = sum(rep(nj ./ aj.^2, mj)) / 2
        Iee = sum(rep((nj .- 1) / σ^2 .+ aj.^-2, mj)) / 2
        Fisher = [Ivv Ive; Ive Iee]
        Vbar = inv(Fisher)

        ### assemble diagonal matrix parts
        g1 = γ .* σ^2 ./ nj
        g3 = zeros(length(nj))
        for l in 1:length(nj)
            g3[l] += sum(Diagonal(nj[l] / aj[l]^4 .*
            [τ^2 * σ^4 .* nj[l] .+ σ^6  -(τ^4 * σ^2 .* nj[l] .+ τ^2 * σ^4);
            -(τ^4 * σ^2 .* nj[l] .+ τ^2 * σ^4) τ^6 .* nj[l] .+ τ^4 * σ^2] * Vbar))
        end

        ### assemble inverse of covariance matrix with Woodbury
        Σdiag = convert(Array{Float64,1}, rep(g1 + 2 * g3, mj))
        Σi = (I - (d ./ Σdiag) * inv(inv(vcov1) + d' * (d ./ Σdiag)) * d') * Diagonal(Σdiag.^-1)

        ###     CONDITIONAL

        ### assemble diagonal matrix parts
        l1 = γ.^2 ./ nj .* σ^2

        ### helpful new parameters
        b = γ ./ nj
        Db = [σ^2; τ^2] / aj.^2

        bvv = - 2 .* nj .* σ^2 ./ aj.^3
        bve = (nj .* τ^2 .- σ^2) ./ aj.^3
        bee = 2 .* τ^2 ./ aj.^3
        D2b = Vector{Vector{Float64}}[ [bvv, bve], [bve, bee]]

        DvIvv = -sum(rep(nj.^3 ./ aj.^3, mj))
        DvIve = sum(rep(nj.^2 ./ aj.^3, mj))
        DvIee = sum(rep(nj / aj.^3, mj))
        DeIvv = -sum(rep(nj.^2 ./ aj.^3, mj))
        DeIve = sum(rep(nj ./ aj.^3, mj))
        DeIee = sum(rep((nj .- 1) / σ^6 + aj.^-3, mj))
        DFisher = Matrix{Float64}[[DvIvv DvIve; DvIve DvIee], [DeIvv DeIve; DeIve DeIee]]

        l3_1 = zeros(length(nj))
        mid = [sum(γ.^2 .* mj) sum(γ.^2 ./ nj .* mj) ] / τ^4
        for l in 1:length(nj)
            for e in 1:2
                l3_1[l] += mid[e] * σ^4 * b[l] * Vbar[e,:]' * Db[:,l]
            end
        end

        l3_2 = zeros(length(nj))
        for l in 1:length(nj)
            for e in 1:2
                for d in 1:2
                    for f in 1:2
                        for g in 1:2
                            l3_2[l] += Vbar[e,f] * Vbar[f,g] * Fisher[e,d] *  b[l] * D2b[e][d][l]
                        end
                    end
                end
            end
        end

        l3_3 = zeros(length(nj))
        for l in 1:length(nj)
            for e in 1:2
                for d in 1:2
                    for f in 1:2
                        l3_3[l] -= 2 * nj[l] * σ^2 * Vbar[e,f] * Fisher[e,d] * b[l] * Vbar[d,:]' * DFisher[e] * Vbar * Db[:,l]
                    end
                end
            end
        end

        l3_4 = zeros(length(nj))
        for l in 1:length(nj)
            for e in 1:2
                for d in 1:2
                    for g in 1:2
                        l3_4[l] += 2 * nj[l] * σ^2 * DFisher[g][e,d] * Vbar[e,d] * b[l] * Vbar[e,:]' * Db[:,l]
                    end
                end
            end
        end
        l3 = l3_1 .+ l3_2 .+ l3_3 .+ l3_4

        l4 = zeros(length(nj))
        for l in 1:length(nj)
            l4[l] += sum(Diagonal(nj[l] / aj[l]^4 .*
            [σ^6  -τ^2 * σ^4; -τ^2 * σ^4 τ^4 * σ^2] * Vbar))
        end

        l5 = zeros(length(nj))
        for l in 1:length(nj)
            DLvv = σ^4
            DLve = - τ^2 * σ^2
            DLee = τ^4
            l5[l] += nj[l] * (σ^2 - 2 * nj[l] * τ^2) * (nj[l] * τ^2 + σ^2)^-4 * sum(Diagonal([DLvv DLve; DLve DLee] * Vbar)) / 2
        end

        ### assemble inverse of covariance matrix with Woodbury
        cΣdiag = convert(Array{Float64,1}, rep(max1(l1 .+ l3 .+ l4 .- l5, 1/100000), mj))
        cΣ = Diagonal(cΣdiag) + d * (vcov1 + vcov1 * [1 sum(xbar); sum(xbar) sum(xbar.^2)] *
                    vcov1 .* sum(nj .* γ .* (1 .- γ) .* mj) ./ σ1^2) * d' + (rep(γ, mj) .* d) *
                    vcov1 * d' + ((rep(γ, mj) .* d) * vcov1 * d')'
        cΣi = inv(cΣ)

        ###     COMPUTE MEAN
        β1 = beta(vcov0, x, xbar, y, ybar, σ1, τ1, nj, mj, γ0)
        μ1 = rep(1 .- γ0, mj) .* (β1[1] .+ β1[2] .* xbar) + rep(γ0, mj) .* ybar

        β2 = beta(vcov1, x, xbar, y, ybar, σ, τ, nj, mj,  γ) # = fixef(model)
        μ2 = rep(1 .- γ, mj) .* (β2[1] .+ β2[2] .* xbar) + rep(γ, mj) .* ybar

        ###     NONCENTRALITY PARAMETER
        ay = d0 * β1 .- rep(1 .- γ0, mj) .* ybar
        λ[i, 1] = maximum([(ay' * cΣ0i * ay)[1] - tr(rep(γ0.^2 ./ nj.^3, mj) .* Diagonal(cΣ0i)) / τ1^4 * σ1^6, 0])

        ay = d * β1 .- rep(1 .- γ, mj) .* ybar
        λ[i, 2] = maximum([(ay' * cΣi * ay)[1] - tr(rep(γ.^2 ./ nj.^3 ./ aj.^2, mj) .* Diagonal(cΣi)) * σ^6 , 0])

        ###     QUANTILES
        qc1 = cquantile(NoncentralChisq(m, λ[i, 1]), 0.05) # quantile of alpha = 0.05
        qc2 = cquantile(NoncentralChisq(m, λ[i, 2]), 0.05) # quantile of alpha = 0.05

        ###     COMPUTE STATISTICS
        t1 = (μ1 - μ0) .+ Δ'
        t2 = (μ2 - μ0) .+ Δ'
        Coverage[1,:] += (sum(t1' .* (t1' * Σ0i), dims = 2) .< qm) / nsim
        Coverage[2,:] += (sum(t2' .* (t2' * Σi), dims = 2) .< qm) / nsim
        tmp1 = sum(t1' .* (t1' * cΣ0i), dims = 2)
        Coverage[3,:] += (tmp1 .< qc0) / nsim
        Coverage[4,:] += (tmp1 .< qc1) / nsim
        Coverage[5,:] += (sum(t2' .* (t2' * cΣi), dims = 2) .< qc2) / nsim

        ###     COMPUTE RELATIVE VOLUME
        # (via approximations as the determinant of large matrices is numerically unstable)
        Volume[i, 3] = prod(Σ0diag ./ cΣ0diag)^-0.5 * (qm / qc0)^(-m / 2)
        Volume[i, 4] = prod(Σ0diag ./ cΣ0diag)^-0.5 * (qm / qc1)^(-m / 2)
        Volume[i, 5] = prod(Σdiag ./ cΣdiag)^-0.5 * (qm / qc2)^(-m / 2)
    end

    ###     SAVE RESULTS
    filename = "data/unbalanced/" * string(now())
    mkpath(filename)
    Specs = DataFrame(nsim = nsim, m = m, ni = nj, mj = mj, n = N, lambda = λ0,
                        sv = τ0, se = σ0, sCv = τ1, sCe = σ1)
    CSV.write(filename * "/Specs.csv",  Specs, writeheader = true)
    CSV.write(filename * "/Coverage.csv",  DataFrame(Coverage), writeheader = false)
    CSV.write(filename * "/Lambda.csv",  DataFrame(λ), writeheader = false)
    CSV.write(filename * "/Volume.csv",  DataFrame(Volume), writeheader = false)

end
