using JLD2
using Interp1d
using Flux
using Flux: crossentropy, onecold, onehotbatch, params, train!
using CSV
using DataFrames
using MAT
include("train_NN.jl")
include("utils_GTM.jl")
include("GTMSSA.jl")
include("constants.jl")
all_flag = true

@load joinpath(MODELWEIGHT_DIR, "model_stats_prob.jld2") model
@load joinpath(DATA_DIR, "test_GTM.jld2") params_arr_test sim_results_pro_test theo_stats_test

params_arr_test = convert(Vector{Vector{Float32}},params_arr_test)
sim_results_pro_test = convert(Vector{Vector{Float32}},sim_results_pro_test)

param_value_arr = params_arr_test
theo_stats = reduce(vcat,transpose.(theo_stats_test))

m_NN = mean.(Distribution.(Ref(model), param_value_arr))
var_NN = var.(Distribution.(Ref(model), param_value_arr))
fano_NN = var_NN./m_NN

mnb = Distribution.(Ref(model), param_value_arr)
skew_NN = StatsBase.skewness.(mnb).+1

mnb = Distribution.(Ref(model), param_value_arr)
kurt_NN = StatsBase.kurtosis.(mnb).+3

mnb = Distribution.(Ref(model), param_value_arr)
bimodcoeff_NN = bimodcoeff_test_.(mnb)

stat_NN = [m_NN fano_NN skew_NN kurt_NN bimodcoeff_NN]
theo_stats_sub = theo_stats[:,[1,3,4,5,6]]
theo_NN_stats = [stat_NN theo_stats_sub]

theo_NN_stats_df = DataFrame(theo_NN_stats, ["m_NN","fano_NN","skew_NN","kurt_NN", "bimodcoeff_NN","m_true","fano_true","skew_true","kurt_true", "bimodcoeff_true"])
CSV.write( joinpath(RESULT_DIR, "theo_NN_stats.csv"),theo_NN_stats_df)

# The following code is to calculate the distribution
function save_prob_result(data,file_path)
    df = DataFrame(data,:auto)
    CSV.write(file_path,df)
end

if all_flag
    X_test = params_arr_test
    y_test = sim_results_pro_test
    # for i=1:length(X_test)
    insert_result_array = Vector{Float64}[]
    for i=1:length(X_test)
        # for i=1:10
        predicted_result_array =  Vector{Float64}[]
        simPra = X_test[i]
        simPra = convert(Vector{Float64},simPra)
        push!(simPra,1)
        ps = simPra[1:5]
        result = gillespieGTM(simPra,i)
        # insert_result = interp_data(result.data[:,3],result.time)
        # proResult = convertCountToProb(insert_result)
        h_reulst = fit(Histogram, result.data[:,3][1:length(result.time)-1],weights(diff(result.time)),nbins=findmax(result.data[:,3])[1])
        proResult = h_reulst.weights./sum(h_reulst.weights)
        predicted_result = pred_pdf(model, ps, 0:length(proResult))

        push!(insert_result_array,insert_result)
        push!(predicted_result_array,predicted_result)
        save_prob_result(predicted_result_array,joinpath(RESULT_DIR, "distributionResult/NN_predicted_$i.csv"))
    end
    save_prob_result(insert_result_array,joinpath(RESULT_DIR, "distributionResult/ssa_result.csv"))
end

@load joinpath("/public/home/zhjiajun/hzw/academic_code/DeepGTM/data", "test_Toy.jld2") params_arr_test sim_results_pro_test theo_stats_test
params_arr_test = convert(Vector{Vector{Float32}},params_arr_test)
sim_results_pro_test = convert(Vector{Vector{Float32}},sim_results_pro_test)
param_value_arr = params_arr_test
X_test = params_arr_test
y_test = sim_results_pro_test
insert_result_array = Vector{Float64}[]
valid_param_list = [10,21,72,106,114,291]
for i in valid_param_list
    predicted_result_array =  Vector{Float64}[]
    simPra = X_test[i]
    simPra = convert(Vector{Float64},simPra)
    push!(simPra,1)
    ps = simPra[1:5]
    result = gillespieGTM(simPra,i)
    # insert_result = interp_data(result.data[:,3],result.time)
    # proResult = convertCountToProb(insert_result)
    h_reulst = fit(Histogram, result.data[:,3][1:length(result.time)-1],weights(diff(result.time)),nbins=findmax(result.data[:,3])[1])
    proResult = h_reulst.weights./sum(h_reulst.weights)
    predicted_result = pred_pdf(model, ps, 0:length(proResult))
    push!(insert_result_array,insert_result)
    push!(predicted_result_array,predicted_result)
    save_prob_result(predicted_result_array,joinpath(RESULT_DIR, "distributionResultFig/NN_predicted_$i.csv"))
end
    save_prob_result(insert_result_array,joinpath(RESULT_DIR, "distributionResultFig/ssa_result.csv"))

i = 2
simPra = X_test[i]
simPra = convert(Vector{Float64},simPra)
push!(simPra,1)
ps = simPra[1:5]
result = gillespieGTM(simPra,i)
insert_result = interp_data(result.data[:,3],result.time)
proResult = convertCountToProb(insert_result)
h_reulst = fit(Histogram, result.data[:,3][1:length(result.time)-1],weights(diff(result.time)),nbins=findmax(result.data[:,3])[1])
h_weight = h_reulst.weights./sum(h_reulst.weights)
normalize(h_reulst)
using Plots
histogram(result.data[1:end-1,3] .- 0.5, bins=-0.5:1:findmax(result.data[:,3])[1], weights=diff(result.time),norm=true,
          size=(600*0.9,400*0.9), grid=false, ticks=true, fontfamily="Helvetica",
          framestyle=:box, xlim=(0,500), legend=false)

histogram(collect(0:236) .- 0.5, bins=-0.5:1:findmax(result.data[:,3])[1], weights=h_weight,norm=true,
          size=(600*0.9,400*0.9), grid=false, ticks=true, fontfamily="Helvetica",
          framestyle=:box, xlim=(0,500), legend=false)

histogram(collect(0:235) .- 0.5, bins=-0.5:1:findmax(result.data[:,3])[1], weights=proResult,norm=true,
size=(600*0.9,400*0.9), grid=false, ticks=true, fontfamily="Helvetica",
framestyle=:box, xlim=(0,500), legend=false)
