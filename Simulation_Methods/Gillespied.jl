# -*- coding: utf-8 -*-
import Pkg
Pkg.add("Random")
Pkg.add("Distribution")
Pkg.add("DelimitedFiles")

using Random, Distribution, DelimitedFiles

function hazzy(x::Vector{Float64}, th::Vector{Float64})::Vector{Float64}
    return [x[1], x[2], x[1], x[2], x[1]].*th[1:5]
end

const post = [[2,0,0,0,1] [0,2,0,0,1]]
const pre = [[1,0,1,0,1] [0,1,0,1,0]]
const S = post - pre
const k = [3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0.0, 8.99e-9, 2e-3]

function gillespied(inits::Vector{Float64}, k::Vector{Float64}, S::Matrix{Int64}, Tmax::Int64, δt::Int64)
    x = inits
    tt = 0.0
    n = trunc(Int, Tmax/δt)
    xmat = fill(-1, (2, n))
    i = 1
    target = 0.0
    C0 = sum(x)
    while i <= n
        error = C0 - sum(x)
        h = hazzy(x, k)
        h0 = sum(h)
        if h0<1e-10
            xmat[:,i:n] = fill(0.0, (2,n-i+1))
            return xmat'
        else
            Exp = Exponential(1/h0)
            tt = tt + rand(Exp)
        end
        while tt>=target && i<=n
            xmat[:,i] = x
            i += 1
            target += δt
        end
        Cat = Categorical(h/h0)
        r = rand(Cat)
        x += S'[:,r]
    end
    return xmat'
end

function replace_nan(x)
    for i=eachindex(x)
        x[i] = isnan(x[i]) ? 0.0 : x[i]
    end
end

function raw_to_summ(sims)::Array{Float64}
    """
    converts the species populations from the gillespie algorithm to 
    copy number and mutation load
    """
    Nsim = size(sims)[3] # no. of simulations
    n = size(sims)[1] # length of one simulation
    out = Array{Float64}(undef, n,2,Nsim)
    for i=1:Nsim
        copy_num = sum(sims[:,:,i], dims=2)
        mut_load = sims[:,2,i]./copy_num
        replace_nan(copy_num)
        replace_nan(mut_load)
        out[:,:,i] = hcat(copy_num, mut_load)
    end
    out
end

function quantiles(sims, p)
    """
    returns quantile summaries from simulations
    """
    Nsim = size(sims)[3] # Nsim: number of simulations
    n = size(sims)[1] # length of one simulation
    out = Array{Float64}(undef, n,length(p),2)
    for t=1:n
        out[t,:,1] = quantile([sims[t,1,i] for i=1:Nsim], p)
        out[t,:,2] = quantile([sims[t,2,i] for i=1:Nsim], p)
    end
    out
end

Nsim = 1000
Tmax = 80*365*24*3600
δt = 24*3600

simulations = Array{Float64}(undef, trunc(Int, Tmax/δt), 2, Nsim)
for i=1:Nsim
    simulations[:,:,i] = gillespied([100.0,100.0], k ,S, Tmax, δt)
end

CH_sim = raw_to_summ(simulations)

quant_sim = quantiles(CH_sim, [0.025,0.1,0.5,0.9,0.975]) 

writedlm("copy_number_quatiles_jl.txt", quant_sim[:,:,1])
writedlm("muation_load_quantiles_jl.txt", quant_sim[:,:,2])












