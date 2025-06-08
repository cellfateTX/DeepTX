struct SSAArgs{X,Ftype,N,P}
    x0::X
    F::Ftype
    nu::N
    parms::P
    tf::Float64
end

struct SSAResult
    time::Vector{Float64}
    data::Matrix{Int64}
    args::SSAArgs
end

function pfsample(w::AbstractArray{Float64,1},s::Float64,n::Int64)
    t = rand() * s
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end

function convertCountToProb(data_counts)
    probArrNum = maximum(data_counts)
    probArrNum = convert(Int64,probArrNum)
    proArr = zeros(probArrNum+1)
    data_counts = convert(Vector{Int64},data_counts)
    for i in eachindex(data_counts)
        # num_index = convert(Int64,data_counts[i]+1)
        num_index = data_counts[i]+1
        proArr[num_index]=proArr[num_index]+1
    end
    proArr=proArr./length(data_counts)
    proArr
end
# function gillespie(x0::AbstractVector{Int64},F::Base.Callable,nu::AbstractMatrix{Int64},parms::Vector{Float64},tf::Float64,i::Int64)
function gillespie(x0::Vector{Int64},F::Base.Callable,nu::Matrix{Int64},parms::Vector{Float64},tf::Float64,i::Int64)

    # Args
    args = SSAArgs(x0,F,nu,parms,tf)
    # Set up time array
    ta = Vector{Float64}()
    t = 0.0
    push!(ta,t)
    # Set up initial x
    nstates = length(x0)
    x = copy(x0')
    xa = copy(Array(x0))
    # Number of propensity functions
    numpf = size(nu,1)
    next_mu = 1
    # Main loop
    # termination_status = "finaltime"
    nsteps = 0
    while t <= tf
        # println("the t is $t")
        if x[1] == 1
            tau_off = rand(Gamma(parms[3],1/parms[4]))
            t_temp = t + tau_off
            while t < t_temp && t < tf
                a0 = parms[6]*x[3];
                if a0==0
                    push!(ta,t)
                    for xx in x + nu[4,:]'
                        push!(xa,xx) 
                    end
                   # update nstep
                    nsteps += 1
                    # x = x + nu[4,:]';
                    break
                end
                tau_1 = rand(Exponential(1/a0))
                # println("the value of tau1 is $tau_1")
                t = t + tau_1
                push!(ta,t)
                x = x + nu[4,:]';
                for xx in x
                    push!(xa,xx)
                end
               # update nstep
               nsteps += 1
            end
            pop!(ta)
            push!(ta,t_temp)
            if t_temp <= tf
                pop!(xa)
                pop!(xa)
                pop!(xa)

                for xx in x - nu[4,:]' + nu[1,:]'
                    push!(xa,xx)
                end
            else
                pop!(xa)
                pop!(xa)
                pop!(xa)
                for xx in x - nu[4,:]'
                    push!(xa,xx)
                end
            end
            x = x - nu[4,:]'+ nu[1,:]' 
            # 把最后时刻重新复制给t。
            t = t_temp
        else
            tau_on = rand(Gamma(parms[1],1/parms[2]))
            t_temp = t + tau_on
            while t < t_temp && t < tf
               amu = [parms[5]*x[2],parms[6]*x[3]]

               a0 = sum(amu) 
               tau_2 = rand(Exponential(1/a0))
               if x[3] == 0
                  next_mu=1
               else
                  next_mu = pfsample(amu,a0,3)
               end
               t = t + tau_2
               push!(ta,t)
               x = x + nu[next_mu+2,:]'
               for xx in x
                   push!(xa,xx)
               end
               nsteps += 1
            end
            pop!(ta)
            push!(ta,t_temp)
            if t_temp <= tf
                pop!(xa)
                pop!(xa)
                pop!(xa)
                if x[3] > 0
                    for xx in x - nu[next_mu+2,:]' + nu[2,:]'
                        push!(xa,xx)
                    end
                else
                    for xx in x - nu[2+2,:]'+ nu[2,:]'
                        push!(xa,xx)
                    end
                end
            else
                pop!(xa) 
                pop!(xa)
                pop!(xa)
                if x[3] > 0
                    for xx in x - nu[next_mu+2,:]'+ nu[2,:]'
                        push!(xa,xx)
                    end
                else
                    for xx in x - nu[2+2,:]'+ nu[2,:]'                        
                        push!(xa,xx)
                    end
                end
            end
            x = x- nu[next_mu+2,:]' + nu[2,:]' 
            t = t_temp
        end
    end
    xar = transpose(reshape(xa,length(x),nsteps+1))
    return SSAResult(ta,xar,args)
end

function F_dd(x,parms)
    (OFF,ON,mRNA) = x
    (kon,koff,ksyn,kdeg) = parms
    [kon*OFF,koff*ON,ksyn*ON,kdeg*mRNA]
end  

function gillespieGTM(parms::Vector{Float64},i;tf=30000.0)
    x0 = [0,1,0]
    nu = [[-1 1 0]; [1 -1 0]; [0 0 1]; [0 0 -1]]
    gillespie(x0,F_dd,nu,parms,tf,i)
end


