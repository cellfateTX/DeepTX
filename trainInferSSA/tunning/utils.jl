using MAT
using Interp1d
using DataFrames
function saveMatFile(data,file_path,file_name)
    sim_param_results_df = DataFrame(data,:auto)
    file = matopen(file_path, "w")
    write(file, file_name, Matrix(sim_param_results_df))
    close(file)
end


function loadMatFile(file_path,var_name)
    
    file = matopen(file_path)
    var_names = keys(file)
    stat_the = read(file, var_name) # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    stat_the
end

function interp_data(yi::Vector{Int64},xi::Vector{Float64})
    yi = convert(Vector{Float64},yi)  
    x = collect(20000:0.05:30000)
    a=Vector{Float32}()
    for mode in INTERP_MODE_LIST
        if string(mode)=="Previous"
            b_index = 1
            f = interp(xi[b_index:length(xi)], yi[b_index:length(xi)], mode) # get an interpolation function
            y = f.(x) # Do interpolation
            a = y
        end
    end
    a
end

function calcuPro(data_counts)
    probArrNum = maximum(data_counts)
    probArrNum = convert(Int64,probArrNum)
    proArr = zeros(probArrNum+1)
    data_counts = convert(Vector{Int64},data_counts)
    for i in eachindex(data_counts)
        # num_index = convert(Int64,data_counts[i]+1)
       
        if data_counts[i]>=0
            num_index = data_counts[i]+1
            proArr[num_index]=proArr[num_index]+1
        end

    end
    proArr=proArr./length(data_counts)
    proArr
    
end

# transform the counts data array to probArr
function convertCountsToProb(data_counts)
    probArrNum = maximum(data_counts)
    probArrNum = convert(Int64, probArrNum)
    proArr = zeros(probArrNum + 1)
    data_counts = convert(Vector{Int64}, data_counts)
    for i in eachindex(data_counts)
        if data_counts[i] >= 0
            num_index = data_counts[i] + 1
            proArr[num_index] = proArr[num_index] + 1
        end
    end
    proArr = proArr ./ length(data_counts)
    proArr
end

# function linear_regression(X::AbstractVector,Y::AbstractVector;log_form=false)
#     if log_form
#         X=log.(X)
#         Y=log.(Y)
#     end
#     data = DataFrame(X=X, Y=Y)
# line = lm(@formula(Y ~ X), data)
# round(r2(line),digits=2) 
# end

function save_prob_result(data,file_path)
    df = DataFrame(data,:auto)
    CSV.write(file_path,df)
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

function slicematrix(A::AbstractMatrix)
    return [A[i, :] for i in 1:size(A,1)]
end

function hellinger2(aa, bb)
    ret = 0.0
    for (a, b) in zip(aa, bb)
        ret += (sqrt(a) - sqrt(b))^2
    end
    ret
end

function calculate_busrt_dynamics(data::DataFrame)
    data[:, :tauOn] = data[!,1]./data[!,2]
    data[:, :tauOff] = data[!,3]./data[!,4]
    data[:, :tauOnRatio] = data[:, :tauOn] ./(data[:, :tauOn]+data[:, :tauOff])
    data[:, :bf] = 1 ./ (data[:, :tauOn] .+ data[:, :tauOff])
    data[:,:bs] = data[:, :tauOn] .* data[!,5]
    data
end



# function for simulation
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

function gillespie(x0::AbstractVector{Int64},F::Base.Callable,nu::AbstractMatrix{Int64},parms::Vector{Float64},tf::Float64,i::Int64)
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
    termination_status = "finaltime"
    nsteps = 0
    while t <= tf
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
            # 新增，要把最后时刻重新复制给t。
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
            # 新增，要把最后时刻重新复制给t。
            t = t_temp
        end
    end
    xar = transpose(reshape(xa,length(x),nsteps+1))
    # print(ta)
    return SSAResult(ta,xar,args)
end

function F_dd(x,parms)
    (OFF,ON,mRNA) = x
    (kon,koff,ksyn,kdeg) = parms
    [kon*OFF,koff*ON,ksyn*ON,kdeg*mRNA]
end  

function insert_data(yi,xi)
    Y = Vector{Vector{Float32}}()
    yi = convert(Vector{Float64},yi)  
    x = collect(1000:0.03:3000)
    a=Vector{Float32}()
    # print(x)
    for mode in INTERP_MODE_LIST
        if string(mode)=="Previous"
            # b_index = convert(Int32,ceil(length(xi)/2))
            b_index = 1
            f = interp(xi[b_index:length(xi)], yi[b_index:length(xi)], mode) # get an interpolation function
            y = f.(x) # Do interpolation
            # p1 = Plots.histogram(y)
            # Plots.savefig("test/fig/myfig$i.png")
            # push!(X, parms[1:(end-1)])
            # push!(Y, y)
            a = y
        end
    end
    # println(a)
    a
end


function calcuPro(data_counts)
    probArrNum = maximum(data_counts)
    probArrNum = convert(Int64,probArrNum)
    proArr = zeros(probArrNum+1)
    data_counts = convert(Vector{Int64},data_counts)
    for i in eachindex(data_counts)
        # num_index = convert(Int64,data_counts[i]+1)
       
        if data_counts[i]>=0
            num_index = data_counts[i]+1
            proArr[num_index]=proArr[num_index]+1
        end

    end
    proArr=proArr./length(data_counts)
    proArr
    
end

function slicematrix(A::AbstractMatrix)
    return [A[i, :] for i in 1:size(A,1)]
end
