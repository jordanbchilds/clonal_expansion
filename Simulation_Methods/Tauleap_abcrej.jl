#!/usr/local/bin/julia
# -*- coding: utf-8 -*-

using Pkg
# Pkg.add("Distributed")
# Pkg.add("Plots")
# Pkg.add("DelimitedFiles")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Random")
# Pkg.add("Distributions")

using Distributed, Plots, DelimitedFiles, CSV, DataFrames, Random, Distributions

length(Sys.cpu_info())
addprocs(24) ;

@everywhere using Random, Distributions

moraes = DataFrame(CSV.File("../moraes_data.csv"))
moraes_df = filter(r->r.time>=15, moraes)

wild_df = filter(row->row.type=="wild", moraes_df)  ; 

wild_mat = Array{Float64, 2}(undef, length(unique(moraes_df.time)), length(unique(moraes_df.bound)))
bounds = unique(wild_df.bound)
for i=1:length(bounds)
        wild_mat[:,i] = filter(r->r.bound==bounds[i], wild_df).exp_level
end

@everywhere struct SPN
    init_pop::Real
    pops::Real
    rate_vec::Vector{Real}
    Stoichiometry_matrix::Vector{Real}
    function SPN(init_pop, pops, rate_vec, Stoichiometry_matrix)
        new(init_pop, pops, rate_vec, Stoichiometry_matrix)
    end
end

@everywhere init(N::SPN) = Float64.(N.init_pop)
@everywhere starting_pop(N::SPN) = Float64.(N.pops)
@everywhere rates(N::SPN) = Float64.(N.rate_vec)
@everywhere StoiMat(N::SPN) = Float64.(N.Stoichiometry_matrix)

@everywhere const post = [2,0] 
@everywhere const pre = [1,1]
@everywhere const S = post - pre

hour = 3600
day = 24*hour
year = 365*day

step_str = "1"
step = 1*day
step_out = 1*day

Tmax = 50*year 
Nsim = 1000 
timed = unique(wild_df[:,"time"])*day ;

@everywhere function hazard(x::Float64, th::Vector{Float64}, error::Float64)::Vector{Float64}
    k = th[1:2]
    Kc = th[3:4]
    if error>=0
        k1 = k[1]+error*Kc[1]
        return x .*[k1, k[2]]
    else 
        k1 = 2*k[1]/(1+exp(-error*Kc[2]))
        return x .*[k1, k[2]]
    end
end

@everywhere function randPois(λ::Vector{Float64})::Vector{Int64}
    pos_rates = (λ.>0).*λ
    [rand(Poisson(rate)) for rate in pos_rates]
end

@everywhere function transform_summ(popdym, C0)::Array{Union{Float64, Missing}}
    popdym / C0
end

@everywhere function tauleap(spn::SPN, Tmax::Real, dt::Real, dtout::Real, target)::Array{Union{Float64, Missing}}

    x = starting_pop(spn) # population to start simulation
    C0 = isnan(init(spn)) ? x : init(spn) # if NaN given for initial population then x starts simuatlion
                                          # if not the given C0 starts the simulation
    k = rates(spn)
    S = StoiMat(spn)
    N = trunc(Int, Tmax/dt) 
    Nout = length(target) 
    popdym = Array{Float64, 2}(undef, length(x),Nout)
    # target = 0.0
    tt = 0.0
    i = 1
    for _=1:N
        while i<=Nout && tt>=target[i] 
            popdym[i] = x
            #target += dtout
            i += 1
        end
        error = C0 - x
        h = hazard(x, k, error)
        if( sum(h) < 1e-10 )
            popdym[i:Nout] = zeros(Nout-i+1)
            return transform_summ(popdym, C0)
        end
        R = randPois(h*dt)
        x = x + (S'*R) 
        x = x<0.0 ? 0.0 : x
        if x > 10*C0
            popdym = fill(NaN, Nout)
            return popdym
        end
        tt += dt
    end
    return transform_summ(popdym, C0)
end


@everywhere function par_sim(Nsim::Int64, f, spn::SPN, Tmax::Real, dt::Real, dtout::Real, target)
    np = nworkers()            # Number of processes available.
    output = Vector{Array{Float64}}(undef, Nsim) 
    i = 1
    nextidx() = (idx = i; i += 1; idx) # Function to know which is the next work item.
    @sync begin #@sync: must complete all jobs in block
        for p = 1:np # loops through all processes (workers)
            if p != myid() || np == 1 # first worker used only if all others are busy 
                @async begin # launch several tasks simultaneaously
                    while true
                        idx = nextidx()
                        if idx > Nsim
                            break
                        end
                        output[idx] = remotecall_fetch(f, p, spn, Tmax, dt, dtout, target)
                    end
                end
            end
        end
    end
    output
end

@everywhere function quantiles(sims, p)::Array{Float64, 2}
    """
    returns quantile summaries from simulations
    """
    Nsim = length(sims) # Nsim: number of simulations
    n = length(sims[1]) # length of one simulation
    np = length(p)
    out = Array{Float64}(undef, n, np)
    for t=1:n
        vec = [ sims[i][t] for i=1:Nsim if !isnan(sims[i][t]) ]
        out[t,:] = length(vec)>0 ? quantile(vec, p) : zeros(np)*NaN
    end
    out
end

function prior_draw()::Vector{Float64}
    # some priors taken from JPMorgan
    k1 = rand( Uniform(2e-8, 8e-6) )
    k2 = rand( Uniform(2e-8, 8e-6) )
    kc1 = rand( Uniform(0, 8e-6) )
    kc2 = rand( Uniform(0, 1e-3) )
    return [k1,k2,kc1,kc2]
end ;

function euclidean_dist(sims::Vector{Array{Float64, 2}}, data::Array{Float64, 2})::Vector{Float64}
    Nabc = length(sims)
    output = Vector{Float64}(undef, Nabc)
    for i=1:Nabc
        output[i] = sum(sims[i] - data)^2
    end
    output
end ;

function easy_abc(Nabc, topn, data, inits, Tmax, dt, dtout, target,  Nsim)
    # draw parameters
    abc_sims = Vector{Array{Float64,2}}(undef, Nabc) # use three as length of  s(z)
    theta_post = Vector{typeof(prior_draw())}(undef, Nabc)
    
    i = 1
    
    while i<=Nabc
        # init_pop, pops, rate_vec, Stoichiometry_matrix
        theta_star = prior_draw()
        spn_star = SPN(inits[1], inits[2], theta_star, inits[3])
        sims = par_sim(Nsim, tauleap, spn_star, Tmax, dt, dtout, target)
        qnts = quantiles(sims, [0.025,0.5,0.975])
        if !isnan(sum(qnts))
            abc_sims[i] = qnts
            theta_post[i] = theta_star
            i += 1
        end        
    end
    dist = euclidean_dist(abc_sims, data)
    topind = Vector{Bool}(undef, Nabc)
    topind .= dist .<= sort(dist)[topn]

    theta_post[topind]
end

Nabc = 10000
Nout = 1000
inits = (200.0, 20.0, S)

@time abc_output = easy_abc(Nabc, Nout, wild_mat, inits, Tmax, step, step_out, timed.-15*day, Nsim) 

mkpath("Inference")
writedlm("Inference/tauleap_abcrej.txt", abc_output)


