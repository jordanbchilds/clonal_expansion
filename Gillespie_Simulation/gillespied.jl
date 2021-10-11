using Random, Distrubtions

function hazard(x, k, error, Kc)
	if error>=0
		k1 = k[1]+error*Kc[1]
		k2 = k[2]+error*Kc[1]
		return [x[1],x[2],x[1],x[2],x[1]].*[k1,k2,k[3],k[4],[5]]
	else
		k1 = 2*k[1]/(1+exp(-error/Kc[2])
		k2 = 2*k[2]/(1+exp(-error/Kc[2])
		return [x[1],x[2],x[1],x[2],x[1]].*[k1,k2,k[3],k[4],[5]]
	end
end

post = [[2,0,0,0,1] [0,2,0,0,1]]
pre = [[1,0,1,0,1] [0,1,0,1,0]]
S = post - pre
k = [3.06e-8, 3.06e-8, 3.06e-8, 3.06e-8, 0.0]

function gillespied(k::Vector{Float64}, S::Matrix{Int64}, Tmax::Int64, dt::Int64)
    # draw intial copy number and mutation load
    copy_num = rand(Normal(200, 50), 1)
    mut_load = rand(Beta(10,18.57), 1)
    # define species population vector
    x = Vector{Float64}()
    append!(x, round.(copy_num.*mut_load, digits=0))
    append!(x, round.(copy_num.*(1 .-mut_load), digits=0))
    
    tt = 0
    n = trunc(Int, Tmax/dt)
    xmat = Matrix{Float64}(undef, (n,2)) # store output
    i = 1
    target = 0
    C_0 = sum(x) # initial cipy number
    
    # repeat indefintely until pulled out of loop
    while true
        h = hazard(x,k) # calc hazards
        h0 = sum(h) 
        Exp = Exponential(h0)
        # Darrens trick to escape simulation if both populations are zero
        if h0<1e-10
	    tt = 1e99
        else
	    Exp = Exponential(h0)
 	    tt += rand(Exp)
	end
        while tt >= target
            xmat[i,:] = x
            i += 1
            target = target + dt
            if i>n # if Tmax reached return output
                return trunc.(Int, xmat) 
            end
        end
        Cat = Categorical(h/h0)
        j = rand(Cat) # which reaciton will occur? j
 	x += S[j,:]   
    end
end

sim1 = gillespied(k, S, 365*24*3600, 24*3600)

open("mtDNA.txt", "w") do file
	write(file, sim1)
end


