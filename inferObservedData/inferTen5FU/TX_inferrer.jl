# infer the param value of the scRNA-seq data by neural network.
using Catalyst
using Random
using BlackBoxOptim
using Flux
using CSV
using DataFrames
using JLD2
using Sobol
include("../train_NN.jl")
include("../utils.jl")
include("../constants.jl")

d = 5
intensity = 1
logranges = [  1.0 15.0
                0.1 10.0
                1.0 15.0
                0.01 10.0
                0.1  400.0
                1.0 1.0
             ]

prior = Product(Uniform.(logranges[1:d, 1], logranges[1:d, 2]))

function hellinger2(aa, bb)
    ret = 0.0
    for (a, b) in zip(aa, bb)
        ret += (sqrt(a) - sqrt(b))^2
    end
    ret
end

function loss_hellinger_map(x::AbstractVector, model, hist_yobs, tt = tt)
    bufs = zeros(Threads.nthreads())
    Threads.@threads for i = 1:length(tt)
        ps = x
        # pred = pred_pdf(model, ps, 0:length(hist_yobs))
        pred = pred_pdf_infe(model, ps, 0:length(hist_yobs))
        bufs[Threads.threadid()] += hellinger2(pred, hist_yobs)
    end
    sum(bufs)
end

function inference_parameters(gene_exp_data,model)
    op4 = names(gene_exp_data)
    print(length(op4))
    estimates = Vector{Float64}[]
    for j in eachindex(op4)
        gene_exp = gene_exp_data[:, j]
        hist_yobs = convertCountsToProb(gene_exp)
        tt = [0, 0]
        yobs = []
        Threads.@threads for i in [1]
            opt_result = bboptimize(
                p -> loss_hellinger_map(p, model, hist_yobs, tt);
                SearchRange = [tuple(logranges[i, :]...) for i = 1:d],
                TraceMode = :silent,
            )
            push!(estimates, best_candidate(opt_result))
            println(best_candidate(opt_result))
        end
    end
    estimates
end

function calculate_bs_bf(estimates)
    burst_freq = 1 ./ (estimates[:,1] ./ estimates[:,2] .+ estimates[:,3] ./ estimates[:,4])
    burst_size = estimates[:,5] .* (estimates[:,1] ./ estimates[:,2])
    mean_val = burst_freq.*burst_size
    burst_freq,burst_size,mean_val
end

function save_csv(estimates,file_path,gene_name_arr)
    estimates_df = DataFrame(estimates, :auto)
    estimates_df = DataFrame(Matrix(estimates_df)', :auto);
    estimates_df.gene_name = gene_name_arr
    CSV.write(file_path, estimates_df)
end

function matrix_to_dataframe(estimates,gene_name_arr)
    estimates_df = DataFrame(estimates, :auto)
    estimates_df = DataFrame(Matrix(estimates_df)', :auto);
    estimates_df.gene_name = gene_name_arr
    estimates_df
end
gene_exp_idu = DataFrame(
    CSV.File(
        joinpath(DATA_DIR,"Ten5FU/doseTen_norm_filter.csv" ),
    ),
)[!,1:10]
gene_exp_dmso = DataFrame(
    CSV.File(
        joinpath(DATA_DIR,"Ten5FU/doseZero_norm_filter.csv" ),
    ),
)[!,1:10]

# load the trained model
@load joinpath(MODELWEIGHT_DIR, "$MODEL_TYPE.jld2") model

estimated_idu = inference_parameters(gene_exp_idu,model)
estimated_dmso = inference_parameters(gene_exp_dmso,model)

estimated_idu = matrix_to_dataframe(estimated_idu,names(gene_exp_idu))
estimated_dmso = matrix_to_dataframe(estimated_dmso,names(gene_exp_idu))

number_col=5
mean_idu = mean.(eachcol(gene_exp_idu))
var_idu = var.(eachcol(gene_exp_idu))

mean_dmso = mean.(eachcol(gene_exp_dmso))
var_dmso = var.(eachcol(gene_exp_dmso))

burst_freq_idu,burst_size_idu,mean_es_idu = calculate_bs_bf(estimated_idu)
estimated_idu.bf=burst_freq_idu
estimated_idu.bs=burst_size_idu
estimated_idu.mean_es = mean_es_idu
estimated_idu.mean_true = mean_idu
estimated_idu.var_true = var_idu

burst_freq_dmso,burst_size_dmso,mean_es_dmso = calculate_bs_bf(estimated_dmso)
estimated_dmso.bf=burst_freq_dmso
estimated_dmso.bs=burst_size_dmso
estimated_dmso.mean_es = mean_es_dmso
estimated_dmso.mean_true = mean_dmso
estimated_dmso.var_true = var_dmso

es_file_name_idu = "Ten5FU/doseTen_estimated_$MODEL_TYPE.csv"
es_file_name_dmso = "Ten5FU/doseZero_estimated_$MODEL_TYPE.csv"
CSV.write(joinpath(RESULT_DIR,"$es_file_name_idu" ), estimated_idu)
CSV.write(joinpath(RESULT_DIR,"$es_file_name_dmso" ), estimated_dmso)
