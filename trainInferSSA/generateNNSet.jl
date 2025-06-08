using Distributions
using DataFrames
using StaticArrays
using Sobol
using JLD2
using StatsBase
using CSV
using MAT
include("utils_GTM.jl")
include("constants.jl")

sim_type = "GTM"
data_path = joinpath(root_path, "data/simulationData")
train_rate = 0.8
valid_rate = 0.1
test_rate = 0.1

@load joinpath(data_path, "params_prob_$sim_type.jld2") sim_param_results sim_prob_results 

stat_the = loadMatFile(joinpath(data_path,"stats_params_$sim_type.mat"),"stat_the")
params_the =  loadMatFile(joinpath(data_path,"stats_params_$sim_type.mat"),"params")
stat_value_arr = slicematrix(Matrix{Float32}(stat_the))

data_length = length(sim_param_results)
train_end_index = floor(Int,data_length*train_rate)
valid_end_index = floor(Int,data_length*(train_rate+valid_rate))

params_arr_train = sim_param_results[1:train_end_index]
sim_results_pro_train = sim_prob_results[1:train_end_index]
theo_stats_train = stat_value_arr[1:train_end_index]
@save joinpath(data_path, "train_$sim_type.jld2") params_arr_train sim_results_pro_train theo_stats_train

params_arr_valid =  sim_param_results[train_end_index+1:valid_end_index]
sim_results_pro_valid = sim_prob_results[train_end_index+1:valid_end_index]
theo_stats_valid = stat_value_arr[train_end_index+1:valid_end_index]
@save joinpath(data_path, "valid_$sim_type.jld2") params_arr_valid sim_results_pro_valid theo_stats_valid

params_arr_test = sim_param_results[valid_end_index+1:end]
sim_results_pro_test = sim_prob_results[valid_end_index+1:end]
theo_stats_test = stat_value_arr[valid_end_index+1:end]
@save joinpath(data_path, "test_$sim_type.jld2") params_arr_test sim_results_pro_test theo_stats_test
