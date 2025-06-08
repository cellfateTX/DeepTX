using Catalyst
using Random
using BlackBoxOptim
using Flux
using CSV
using DataFrames
using JLD2
using Sobol
include("train_NN.jl")
include("utils_GTM.jl")

model_type = "GTM"
d = 5
intensity = 1
logranges = [  1.0 15.0
                0.1 10.0
                1.0 15.0
                0.01 10.0
                0.1  400.0
             ]

prior = Product(Uniform.(logranges[1:d, 1], logranges[1:d, 2]))

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
function loss(param, hist_yobs,model)
    pred = pred_pdf_infe(model, param, 0:length(hist_yobs))
    loss =  hellinger2(pred, hist_yobs)
    return hellinger2(pred, hist_yobs)
end
function TX_inferrer(hist_yobs,model;param=[], η=0.01, 
    max_epoch=50, patience=5, min_delta=1e-4)

    logranges = [
        1.0 15.0;
        0.1 10.0;
        1.0 15.0;
        0.01 10.0;
        0.1 400.0
    ]
    
    best_loss = Inf
    wait = 0
    temp_loss = Inf
    if length(param)==0
        prior = Product(Uniform.(logranges[1:d, 1], logranges[1:d, 2]))
        param = rand(prior)
    end
    for epoch in 1:max_epoch
        grads = gradient(p -> loss(p, hist_yobs, model), param)
        grad_values = grads[1]
        new_param = param .- η .* grad_values
        for i in 1:length(param)
            new_param[i] = clamp(new_param[i], logranges[i,1], logranges[i,2])
        end
        param .= new_param
        pred = pred_pdf_infe(model, param, 0:length(hist_yobs))
        temp_loss = hellinger2(pred, hist_yobs)
    end
    pred = pred_pdf_infe(model, param, 0:length(hist_yobs))
    temp_loss = hellinger2(pred, hist_yobs)
    return param, temp_loss
end

# infer the param value of the test set by neural network.
@load joinpath(DATA_DIR, "test_$model_type.jld2") params_arr_test sim_results_pro_test
@load joinpath(MODELWEIGHT_DIR, "model_stats_prob.jld2") model

estimates_BB = []
s = SobolSeq(logranges[:,1] ,logranges[:,2])
ps_train = [ Sobol.next!(s) for i in 1: 100 ]
start_point_num = 1
params_num = length(sim_results_pro_test)
# params_num = 10

# The BBoptimize method provides a good initial value for subsequent gradient descent optimization. 
for i in eachindex(sim_results_pro_test[1:params_num])
    hist_yobs = sim_results_pro_test[i]
    tt = [0,0]
    yobs = []
    p = Progress(1)
    Threads.@threads for i in 1:start_point_num  
        opt_result = bboptimize(p -> loss_hellinger_map(p, model,hist_yobs,tt),ps_train[i]; SearchRange = [ tuple(logranges[i,:]...) for i in 1:d ], TraceMode=:silent)
        push!(estimates_BB, best_candidate(opt_result))
        print(best_candidate(opt_result))
        ProgressMeter.next!(p)
    end
end


estimates = []
losses = []

for i in 1:length(estimates_BB)
    hist_yobs = sim_results_pro_test[i]
    blackBox_param = deepcopy(estimates_BB[i])
    param, temp_loss = TX_inferrer(hist_yobs, model; param=blackBox_param)
    push!(estimates, param)
    push!(losses, temp_loss)
end



@save joinpath(RESULT_DIR, "synthetic/estimated_params_test_$model_type.jld2") estimates

burst_freq_es_li = Float64[]
burst_freq_sim_li = Float64[]
burst_size_es_li = Float64[]
burst_size_sim_li = Float64[]

for i in eachindex(estimates)
    # get the burst frequency from inferred params
    burst_freq_es = 1/(estimates[i][1]/estimates[i][2]+estimates[i][3]/estimates[i][4])
    push!(burst_freq_es_li,burst_freq_es)

    # get the burst frequency from simulated params
    burst_freq_sim =  1/(params_arr_test[i][1]/params_arr_test[i][2]+params_arr_test[i][3]/params_arr_test[i][4])
    push!(burst_freq_sim_li,burst_freq_sim)

    # get the burst size from inferred params
    burst_size_es = estimates[i][5]*(estimates[i][1]/estimates[i][2])
    push!(burst_size_es_li,burst_size_es)

    # get the burst size from simulated params
    burst_size_sim =  params_arr_test[i][5]*(params_arr_test[i][1]/params_arr_test[i][2])
    push!(burst_size_sim_li,burst_size_sim)
end

mean_es = burst_freq_es_li.*burst_size_es_li
mean_sim = burst_freq_sim_li.*burst_size_sim_li

estimates_matrix = reduce(vcat,transpose.(estimates))
paramsArr_matrix = reduce(vcat,transpose.(params_arr_test[1:params_num]))

tau_on_es =  estimates_matrix[:,1]./estimates_matrix[:,2]
tau_on_true = paramsArr_matrix[:,1]./paramsArr_matrix[:,2]

tau_off_es =  estimates_matrix[:,3]./estimates_matrix[:,4]
tau_off_true = paramsArr_matrix[:,3]./paramsArr_matrix[:,4]

paramsArr_statics_matrix_true = [paramsArr_matrix tau_on_true tau_off_true burst_freq_sim_li burst_size_sim_li mean_sim]
paramsArr_statics_matrix_es = [estimates_matrix tau_on_es tau_off_es burst_freq_es_li burst_size_es_li mean_es]

paramsArr_statics_true_es = [paramsArr_statics_matrix_true paramsArr_statics_matrix_es]

theo_NN_estimated_param = DataFrame(paramsArr_statics_true_es, 
    ["ron_true","kon_true","roff_true","koff_true", "mu_true","tau_on_true","tau_off_true","bf_true","bs_true","mean_true",
     "ron_es","kon_es","roff_es","koff_es", "mu_es","tau_on_es","tau_off_es","bf_es","bs_es","mean_es"])
CSV.write( joinpath(RESULT_DIR, "synthetic/theo_NN_estimated_param_$model_type.csv"),theo_NN_estimated_param)
