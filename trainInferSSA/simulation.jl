using Distributions
using DataFrames
using StaticArrays
using Sobol
using Interp1d
using JLD2
using StatsBase
using CSV
using MAT
include("GTMSSA.jl")
include("utils_GTM.jl")
include("constants.jl")

# The Mat parameter is used to verify the effectiveness of the algorithm.
data_path = joinpath(root_path, "data/simulationData")

sim_type = "GTM"
sobel_flag = true

sim_size = 20
train_rate = 0.8
valid_rate = 0.1
test_rate = 0.1
# Parameter range of algorithm simulation
ranges = [  1.0 15.0
            0.1 10.0
            1.0 15.0
            0.01 10.0
            0.1  400.0
            1.0 1.0
            ]

s = SobolSeq(ranges[:,1], ranges[:,2])
ps_train = [ Sobol.next!(s) for i in 1: sim_size ]
colun_num = size(ps_train,1)

sim_param_results=Vector{Float64}[]
sim_prob_results=Vector{Float64}[]
sim_count_results=Vector{Float64}[]

for i = 1:colun_num
    params= ps_train[i]
    try
        @time result = gillespieGTM(params,i)
        interp_result = interp_data(result.data[:,3],result.time)
        push!(sim_count_results, interp_result)
        proResult = convertCountToProb(interp_result)
        push!(sim_prob_results, proResult)
        push!(sim_param_results, params[1:5])
    catch e
        print("error")
        print(e)
        continue
    end
end

#Store simulation parameters and use them to calculate theoretical statistical values ​​in matlab code
saveMatFile(sim_param_results,joinpath(DATA_DIR, "params_$sim_type.mat"),"sim_params$sim_type")
@save joinpath(DATA_DIR, "params_prob_$sim_type.jld2") sim_param_results sim_prob_results
