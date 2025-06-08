using JLD2
using MAT
using Catalyst
using DataFrames
using CSV
using Random
using BlackBoxOptim
using Flux
using Sobol

include("train_NN.jl")
include("utils_GTM.jl")
include("constants.jl")
d = 5
intensity = 1
result_dir = joinpath(root_path, "result/synthetic")

logranges = [  1.0 15.0
                0.1 10.0
                1.0 15.0
                0.01 10.0
                0.1  400.0
             ]

@load joinpath(MODELWEIGHT_DIR, "model_stats_prob.jld2") model

function loss_hellinger_map(x::AbstractVector, model,hist_yobs, tt=tt)
    bufs = zeros(Threads.nthreads())
    Threads.@threads for i in 1:length(tt)
        ps = x
        # pred = pred_pdf(model, ps, 0:length(hist_yobs))
        pred = pred_pdf_infe(model, ps, 0:length(hist_yobs))
        bufs[Threads.threadid()] += hellinger2(pred, hist_yobs)
    end
    sum(bufs)
end

function generate_sub_prob(sim_count,sample_size)
    sim_count_results = convert(Vector{Vector{Int}},sim_count)
    sample_index = rand(1:length(sim_count[1]),sample_size)
    sim_count_results_subset = [sim_count_el[sample_index] for sim_count_el in sim_count_results]
    prob_data = convertCountsToProb.(sim_count_results_subset)
    prob_data = convert(Vector{Vector{Float32}},prob_data)
    prob_data
end

sim_type = "tenFold_1"
nsamples = Int.([ 1e2,5e2,1e3,1.5e3,2e3])
@load joinpath(DATA_DIR, "train_$sim_type.jld2") params_arr_train sim_results_pro_train theo_stats_train sim_count_train
params_arr_test = params_arr_train

@load joinpath(DATA_DIR, "test_$sim_type.jld2") params_arr_test sim_results_pro_test theo_stats_test sim_count_test
sim_count_train = sim_count_test

params_num = 2
for nsample in nsamples
    s = SobolSeq(logranges[:,1] ,logranges[:,2])
    ps_train = [ Sobol.next!(s) for i in 1: 100 ]
    start_point_num = 30
    for j in eachindex(sim_count_train[1:params_num])
        estimates = []
        tt = [0,0]
        yobs = []
        p = Progress(1)
        Threads.@threads for i in 1:start_point_num  
            sim_results_pro_test = generate_sub_prob(sim_count_train[1:params_num],nsample)
            hist_yobs = sim_results_pro_test[j]
            opt_result = bboptimize(p -> loss_hellinger_map(p, model,hist_yobs,tt),ps_train[i]; SearchRange = [ tuple(logranges[i,:]...) for i in 1:d ], TraceMode=:silent)
            push!(estimates, best_candidate(opt_result))
            print(best_candidate(opt_result))
            ProgressMeter.next!(p)
        end
    
        estimates_matrix = reduce(vcat,transpose.(estimates))
        # paramsArr_matrix = reduce(vcat,transpose.(params_arr_test[i]))

        tau_on_es =  estimates_matrix[:,1]./estimates_matrix[:,2]
        tau_off_es =  estimates_matrix[:,3]./estimates_matrix[:,4]
        burst_freq_es = 1 ./(tau_on_es .+ tau_off_es)
        burst_size_es = estimates_matrix[:,5] .* tau_on_es
        mean_es = burst_freq_es.*burst_size_es

        # paramsArr_statics_true_es = [paramsArr_statics_matrix_true paramsArr_statics_matrix_es]
        paramsArr_statics_matrix_es = [estimates_matrix tau_on_es tau_off_es burst_freq_es burst_size_es mean_es]
        theo_NN_estimated_param_es = DataFrame(paramsArr_statics_matrix_es, 
        ["ron_es","kon_es","roff_es","koff_es", "mu_es","tau_on_es","tau_off_es","bf_es","bs_es","mean_es"])
        # theo_NN_estimated_param = DataFrame(paramsArr_statics_true_es, 
        #     ["ron_true","kon_true","roff_true","koff_true", "mu_true","tau_on_true","tau_off_true","bf_true","bs_true","mean_true",
        #     "ron_es","kon_es","roff_es","koff_es", "mu_es","tau_on_es","tau_off_es","bf_es","bs_es","mean_es"])
        CSV.write( joinpath(result_dir, "validation/theo_NN_estimated_param_es2_$j-$nsample.csv"),theo_NN_estimated_param_es)

    end
end

paramsArr_matrix = reduce(vcat,transpose.(params_arr_test[1:params_num]))
tau_on_true = paramsArr_matrix[:,1]./paramsArr_matrix[:,2]
tau_off_true = paramsArr_matrix[:,3]./paramsArr_matrix[:,4]
burst_freq_sim =  1 ./(tau_on_true .+ tau_off_true)
burst_size_sim = paramsArr_matrix[:,5] .* tau_on_true
mean_sim = burst_freq_sim.*burst_size_sim
paramsArr_matrix = reduce(vcat,transpose.(params_arr_test[1:params_num]))
paramsArr_statics_matrix_true = [paramsArr_matrix tau_on_true tau_off_true burst_freq_sim burst_size_sim mean_sim]
theo_NN_estimated_param_true = DataFrame(paramsArr_statics_matrix_true, 
["ron_true","kon_true","roff_true","koff_true", "mu_true","tau_on_true","tau_off_true","bf_true","bs_true","mean_true"])
CSV.write( joinpath(result_dir, "validation/theo_NN_estimated_param_true.csv"),theo_NN_estimated_param_true)
