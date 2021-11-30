#!/usr/bin/julia

import Pkg
using Plots, DelimitedFiles, Colors, Random, Statistics
using Distributed

length(Sys.cpu_info())
addprocs(24) ;

@everywhere using Random, Distributed

Tmax = 80*365
Δtout = 7 # time step to save output
Δt = 1 # time step for the simulation
Nsim = 500 ; # number of simulations

@everywhere struct mtDNA 
    rates::Vector{Real}
    unique_id::Int
    parent_id::Int
    status::String
    
    function mtDNA(rates, unique_id, parent_id, status)
        if !(status in ["wild", "mutant"])
            error("Molecules must be of type 'wild' or 'mutant' ")
        end
        if status=="wild"
            new(rates, unique_id, parent_id, status)
        elseif status=="mutant"
            new(rates, unique_id, parent_id, status)
        end
    end
end 

@everywhere rates(mol::mtDNA) = mol.rates
@everywhere unique_id(mol::mtDNA) = mol.unique_id
@everywhere parent_id(mol::mtDNA) = mol.parent_id
@everywhere status(mol::mtDNA) = mol.status 

@everywhere function counter(system_state)::Vector{Int64}
    """
    Calculates the population size for wild and mutant type
    """
    copy_num = length(system_state)
    W = sum([1 for x in system_state if status(x)=="wild"])
    return [W, copy_num-W]
end

@everywhere function transform_summ(pop_dynamics)::Array{Any}
    copy_num = pop_dynamics[:,1] .+ pop_dynamics[:,2]
    mut_load = pop_dynamics[:,2]./copy_num
    return hcat(copy_num, mut_load)
end

@everywhere function agented(init, Tmax::Real, dt::Real, out_dt::Real)
    N = trunc(Int, Tmax/dt)
    Nout = trunc(Int, Tmax/out_dt)
    system_state = init
    current_id = length(init) + 1
    popdym = Array{Union{Float64, Missing}}(undef, (2,Nout+1))
    C0 = length(init)
    target = 0.0
    tt = 0.0
    i = 1
    for k=1:N
        if tt>=target
            popdym[:,i] = counter(system_state)
            target += out_dt
            i += 1
        end
        molecules_to_remove = Vector{Int}()
        new_molecules = Vector{mtDNA}()
        for mol_ind=1:length(system_state)
            molecule = system_state[mol_ind]
            roll = rand(Float64)
            cdf = cumsum( rates(molecule) ) 
            if 0.0<roll && roll<cdf[1] # degredation
                append!( molecules_to_remove, mol_ind)
            elseif cdf[1]<roll && roll<cdf[2] # replication
                append!(molecules_to_remove, mol_ind)
                for j=1:2
                    current_id += 1
                    daughter = mtDNA([3.06e-8,3.06e-8,0]*3600*24, current_id, unique_id(molecule), status(molecule) )
                    push!(new_molecules, daughter)
                end
            elseif cdf[2]<roll && roll<=cdf[3] # mutation
            # mutation last as has smallest probability
               append!(molecules_to_remove, mol_ind)
                for j=1:2
                    current_id += 1
                    daughter = mtDNA([3.06e-8,3.06e-8,0]*3600*24, current_id, unique_id(molecule), ["wild","mutant"][j])
                    push!(new_molecules, daughter)
                end
            end
        end
        system_state = [mol for (i,mol) in enumerate(system_state) if !(i in molecules_to_remove) ]
        append!(system_state, new_molecules)
        tt += dt
        if sum(counter(system_state)) == 0
            popdym[:,i:Nout] = fill(missing, (2,Nout-i+1))
            popdym = popdym'
            copy_num = popdym[:,1] .+ popdym[:,2]
            mut_load = popdym[:,2] ./ copy_num
            return hcat(copy_num, mut_load)
        end
    end
    popdym = popdym'
    copy_num = popdym[:,1] .+ popdym[:,2]
    mut_load = popdym[:,2] ./ copy_num
    return hcat(copy_num, mut_load)
end

C0 = 200
h = 0.5
W0 = round.( C0.*(1 .-h), digits=0)
M0 = round.( C0.*h, digits=0)
inits = [mtDNA([3.06e-8, 3.06e-8, 0]*3600*24, x,-1,"wild") for x=1:W0 ]# initial state of system
append!(inits, [mtDNA([3.06e-8, 3.06e-8, 0]*3600*24, x,-1,"mutant") for x=W0+1:W0+M0] ) ;

@time abm_sim = agented(inits, Tmax, Δt, Δtout)

"""
one simple simulation takes 0 - 3 seconds
"""

abm_sim ;

# @time map(agented, Nlist, Tmaxs, Δts, Δtouts) 
"""
500 simple simulations: 950 seconds 
"""

function par_map(Nsim, f, init, Tmax, Δt, Δtout)
    np = nprocs()  # determine the number of processes available
    results = Vector{Any}(undef, Nsim)
    i = 1
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx > Nsim
                            break
                        end
                        results[idx] = remotecall_fetch(f, p, inits, Tmax, Δt, Δtout)
                    end
                end
            end
        end
    end
    results
end

@time simulations = par_map(Nsim, agented, inits, Tmax, Δt, Δtout) ; 
"""
simple simulution : 280 seconds
4 workers
time step = 1 day

simple simulation : [a long fucking time] 6700 seconds
(legit just 24 times longer than the day)
4 workers
time step = 1 hour
"""

function quantiles(sims, p)
    """
    returns quantile summaries from simulations
    """
    Nsim = length(sims) # Nsim: number of simulations
    n = size(sims[1])[1] # length of one simulation
    out = Array{Float64}(undef, n,length(p),2)
    for t=1:n
        out[t,:,1] = quantile(skipmissing([sims[i][t,1] for i=1:Nsim]), p)
        out[t,:,2] = quantile(skipmissing([sims[i][t,2] for i=1:Nsim]), p)
    end
    out
end

sims_qntl = quantiles(simulations, [0.025,0.25,0.5,0.75,0.975]) ;

myBlack = colorant"rgb(0,0,0,0.1)"
ts = [0:Δtout:Tmax;]./(365);

typeof(simulations)
n = trunc(Int, Tmax/Δtout)
sim_mat = Array{Union{Float64, Missing}}(undef, (n+1,2,Nsim))

for i=1:Nsim
    sim_mat[:,:,i] = simulations[i]
end

p1 = plot(ts, tt[:,1,:], color=myBlack, legend=false, title="Copy Number")
p2 = plot(ts, tt[:,2,:], color=myBlack, legend=false, title="Mutation Load")
plot(p1, p2, layout=(1,2), legend=false)
savefig("Simulations/PDF/abm_simulations_oneday.pdf")

p3 = plot(ts, sims_qntl[:,:,1], title="Copy Number Qunatiles")
p4 = plot(ts, sims_qntl[:,:,2], title="Mutation Load Quantiles")
plot(p3, p4, layout=(1,2), legend=false)
savefig("Simulations/PDF/abm_qntls_oneday.pdf")

writedlm("Simulations/CN_qnt_abm_jl_oneday.txt", sims_qntl[:,:,1])
writedlm("Simulations/ML_qnt_abm_jl_oneday.txt", sims_qntl[:,:,2])


