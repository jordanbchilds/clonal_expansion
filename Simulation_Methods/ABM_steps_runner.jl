#!usr/local/bin/julia

import Pkg
using Plots, DelimitedFiles, Colors, Random, Statistics
using Distributed

length(Sys.cpu_info())
addprocs(24) ;

@everywhere using Random, Distributed, KernelDensity

const hour = 3600
const day = 24*hour
const year = 365*day

#step_str = "50"

#step = 50*day # time step for the simulation
step_out = 7*day # time step to save output
Tmax = 80*year

Nsim = 1000 ; # number of simulations

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

@everywhere function counter(system_state::Vector{mtDNA})::Vector{Float64}
    """
    Calculates the population size for wild and mutant type
    """
    copy_num = length(system_state)
    W = sum([1 for mol in system_state if status(mol)=="wild"])
    return [W, copy_num-W]
end

@everywhere function transform_summ(popdym)::Array{Float64, 2}
    """
    Transforms the population sizes in to copy number and mutation load
    """
    copy_num = popdym[:,1] .+ popdym[:,2]
    
    mut_load = zeros(length(copy_num))*NaN
    nonzero_idx = copy_num .!= 0.0
    mut_load[nonzero_idx] .=  popdym[nonzero_idx,2] ./ copy_num[nonzero_idx]
    
    return hcat(copy_num, mut_load)
end

@everywhere function agented(init::Vector{mtDNA}, Tmax::Real, dt::Real, dtout::Real)::Array{Float64, 2}
    N = trunc(Int, Tmax/dt) + 1
    Nout = trunc(Int, Tmax/dtout) + 1
    system_state = init
    current_id = length(init) + 1
    popdym = Array{Float64}(undef, (2,Nout))
    C0 = length(init)
    
    target = 0.0
    tt = 0.0
    i = 1
    for _=1:N
        while tt>=target && i<=Nout
            popdym[:,i] = counter(system_state)
            target += dtout
            i += 1
        end
        
        molecules_to_remove = Vector{Int}()
        new_molecules = Vector{mtDNA}()
        for mol_ind=1:length(system_state)
            molecule = system_state[mol_ind]
            roll = rand(Float64)
            cdf = cumsum( rates(molecule) ) / sum( rates(molecule) )
            if 0.0<roll && roll<cdf[1] # degredation
                append!( molecules_to_remove, mol_ind )
            elseif cdf[1]<roll && roll<cdf[2] # replication
                append!(molecules_to_remove, mol_ind)
                for j=1:2
                    current_id += 1
                    daughter = mtDNA([3.06e-8,3.06e-8,0], current_id, unique_id(molecule), status(molecule) )
                    push!(new_molecules, daughter)
                end
            elseif cdf[2]<roll && roll<=cdf[3] # mutation
            # mutation last as has smallest probability
               append!(molecules_to_remove, mol_ind)
                for j=1:2
                    current_id += 1
                    daughter = mtDNA([3.06e-8,3.06e-8,0], current_id, unique_id(molecule), ["wild","mutant"][j])
                    push!(new_molecules, daughter)
                end
            end
        end
        system_state = [mol for (i,mol) in enumerate(system_state) if !(i in molecules_to_remove) ]
        append!(system_state, new_molecules)
        tt += dt
        if length(system_state) == 0
            popdym[:,i:Nout] = zeros(2, Nout-i+1)
            return transform_summ(popdym')
        end
    end
    return transform_summ(popdym')
end

C0 = 200
h = 0.5
W0 = round.( C0.*(1 .-h), digits=0)
M0 = round.( C0.*h, digits=0)
inits = [mtDNA([3.06e-8, 3.06e-8, 0], x,-1,"wild") for x=1:W0 ]# initial state of system
append!(inits, [mtDNA([3.06e-8, 3.06e-8, 0], x,-1,"mutant") for x=W0+1:W0+M0] ) ;

# @time map(agented, Nlist, Tmaxs, Δts, Δtouts) 
"""
500 simple simulations: 950 seconds 
"""

function par_map(Nsim, f, init, Tmax, dt, dtout)
    np = nworkers()  # determine the number of processes available
    results = Vector{Array{Union{Float64, Missing}}}(undef, Nsim)
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
                        results[idx] = remotecall_fetch(f, p, inits, Tmax, dt, dtout)
                    end
                end
            end
        end
    end
    results
end

# @time simulations = par_map(Nsim, agented, inits, Tmax, step, step_out) ; 

"""
time step: 50 days
time taken: 106 seconds
4 workers, 1000 simulations
"""

function quantiles(sims, p)
    """
    returns quantile summaries from simulations
    """
    Nsim = length(sims) # Nsim: number of simulations
    n = size(sims[1])[1] # length of one simulation
    out = Array{Float64}(undef, n,length(p),2)
    for t=1:n
        out[t,:,1] = quantile([sims[i][t,1] for i=1:Nsim if !isnan(sims[i][t,1])], p)
        out[t,:,2] = quantile([sims[i][t,2] for i=1:Nsim if !isnan(sims[i][t,2])], p)
    end
    out
end

@everywhere function par_times(Nsim, f, inits, Tmax, step, step_out)
    np = nworkers()            # Number of processes available.
    Nout = trunc(Int, Tmax/step_out) + 1 # dimension for output
    output = Vector{Vector{Union{Float64, Array{Float64,2}}}}(undef, Nsim) # Where we will write the results. As we do not know
                             # the type (Integer, Tuple...) we write "Any"
    i = 1
    nextidx() = (idx = i; i += 1; idx) # Function to know which is the next work item.
                                       # In this case it is just an index.
    @sync begin #@sync: must complete all jobs in block
        for p = 1:np # loops through all processes (workers)
            if p != myid() || np == 1 # first worker used only if all others are busy 
                @async begin # launch several tasks simultaneaously
                    while true
                        idx = nextidx()
                        if idx > Nsim
                            break
                        end
                        output[idx] = Vector{Union{Float64, Array{Float64,2}}}(undef, 2)
                        output[idx][2] = @elapsed sim = remotecall_fetch(f, p, inits, Tmax, step, step_out)
                        output[idx][1] = sim
                    end
                end
            end
        end
    end
    output
end

step_days = [1,2,5,10,15,25,50]
times = zeros(length(step_days))

myBlack = colorant"rgb(0,0,0,0.1)"
ts = [0:step_out:Tmax;]./year;

for  i=1:length(step_days)
    step_str = string(step_days[i])
    step = step_days[i]*day
    
    simulations_times = par_times(Nsim, agented, inits, Tmax, step, step_out)

    simulations = Vector{Array{Float64,2}}(undef, Nsim)
    sim_times = Vector{Float64}(undef, Nsim)
    for i=1:Nsim
        simulations[i] = simulations_times[i][1]
        sim_times[i] = simulations_times[i][2]
    end


    sims_qntl = quantiles(simulations, [0.025,0.25,0.5,0.75,0.975])
    
    mkpath("Simulations")
    writedlm(string("Simulations/CN_qnt_abm_",step_str,"d.txt"), sims_qntl[:,:,1])
    writedlm(string("Simulations/ML_qnt_abm_",step_str,"d.txt"), sims_qntl[:,:,2])
    writedlm(string("Simulations/abm_",step_str,"d_times.txt"), sim_times)

end




