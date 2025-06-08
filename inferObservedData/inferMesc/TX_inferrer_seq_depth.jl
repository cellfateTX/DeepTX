# infer the param value of the scRNA-seq data by neural network.
using Catalyst
using Random
using BlackBoxOptim
using Flux
using CSV
using DataFrames
using JLD2
using Sobol
using LinearAlgebra
include("../train_NN.jl")
include("../utils.jl")
include("../GTMSSA.jl")
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

function loss(param, hist_yobs,model)
    pred = pred_pdf_infe(model, param, 0:length(hist_yobs))
    return hellinger2(pred, hist_yobs)
end

function TX_inferrer(hist_yobs,model;param=[], η=0.01, 
    max_epoch=50, patience=5, min_delta=1e-4,sample_id=1)
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
        # grad_values .= map(x -> x > 0 ? clamp(x, 0.005, 10.0) : clamp(x, -10.0, -0.005), grad_values)

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

function TX_inferrer_fine_tune(estimates_BB,gene_exp_data)
    estimates = []
    losses = []
    for i in 1:length(estimates_BB)
        hist_yobs = convertCountsToProb(gene_exp_data[!,i])
        blackBox_param = deepcopy(estimates_BB[i])
        param, temp_loss = TX_inferrer(hist_yobs, model; param=blackBox_param)
        push!(estimates, param)
        push!(losses, temp_loss)
    end 
    estimates 
end

function Hierarchical_model(a,seq_depth)
    sample_result = Vector{Int32}() 
    for i in 1:10000
        mRNA_counts_ran = [0:length(a)-1...]
        result = wsample(mRNA_counts_ran, a,1)
        point = rand(Binomial(result[1], seq_depth),1)
        push!(sample_result,point[1])
    end
    sample_result_prob = calcuPro(sample_result)
    sample_result_prob
end

function hellinger_distance(p::Vector{Float64}, q::Vector{Float64})
    # 对齐长度
    n = max(length(p), length(q))
    p_padded = vcat(p, zeros(n - length(p)))
    q_padded = vcat(q, zeros(n - length(q)))
    
    # 归一化为概率分布
    p_norm = p_padded / sum(p_padded)
    q_norm = q_padded / sum(q_padded)

    # 计算 Hellinger 距离
    return norm(sqrt.(p_norm) .- sqrt.(q_norm)) / sqrt(2)
end

gene_exp_idu = DataFrame(
    CSV.File(
        joinpath(DATA_DIR,"mesc/filtered_mesc_data.csv" ),
    ),
)[!,2:end]
gene_exp_idu = DataFrame(permutedims(Matrix(gene_exp_idu)), :auto)[!,1:100]


# load the trained model
@load joinpath(MODELWEIGHT_DIR, "$MODEL_TYPE.jld2") model

estimated_idu_BB = inference_parameters(gene_exp_idu,model)
estimated_idu = TX_inferrer_fine_tune(estimated_idu_BB,gene_exp_idu)
estimated_idu = matrix_to_dataframe(estimated_idu,names(gene_exp_idu))


mean_idu = mean.(eachcol(gene_exp_idu))
var_idu = var.(eachcol(gene_exp_idu))
burst_freq_idu,burst_size_idu,mean_es_idu = calculate_bs_bf(estimated_idu)
estimated_idu.bf=burst_freq_idu
estimated_idu.bs=burst_size_idu
estimated_idu.mean_es = mean_es_idu
estimated_idu.mean_true = mean_idu
estimated_idu.var_true = var_idu

burst_seqDepth_file = "mesc/burst_$SEQ_DEPTH.csv"
CSV.write(joinpath(RESULT_DIR,"$burst_seqDepth_file" ), estimated_idu)

burst_seqDepth_file = "mesc/burst_0.5.csv"
estimated_idu2 = CSV.read(joinpath(RESULT_DIR,"$burst_seqDepth_file" ), DataFrame)

burst_seqDepth_file = "mesc/burst_0.3.csv"
estimated_idu3 = CSV.read(joinpath(RESULT_DIR,"$burst_seqDepth_file" ), DataFrame)

cor(estimated_idu2[!,"bs"],estimated_idu3[!,"bs"])



x_list = [estimated_idu2[!,"bs"],estimated_idu2[!,"bf"], estimated_idu2[!,"bs"],estimated_idu2[!,"bf"]]
y_list = [estimated_idu[!,"bs"], estimated_idu[!,"bf"], estimated_idu3[!,"bs"], estimated_idu3[!,"bf"]]
# 创建一个4行1列的布局
using Plots
layout = @layout [a b c d]

# 准备子图集合
plots = []

# 遍历每一对向量，创建散点图并计算相关系数
for i in 1:4
    x = x_list[i]
    y = y_list[i]
    r = cor(x, y)  # 皮尔逊相关系数
    p = scatter(x, y,
        title = "Scatter plot $i\nr = $(round(r, digits=3))",
        xlabel = "x$i",
        ylabel = "y$i",
        legend = false,
        markersize = 3
    )
    push!(plots, p)
end

# 合并子图并显示
plot(plots..., layout = layout, size = (1600, 400))




