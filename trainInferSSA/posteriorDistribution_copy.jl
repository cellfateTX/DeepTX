using Catalyst
using Random
using BlackBoxOptim
using Flux
using CSV
using DataFrames
using JLD2
using Sobol
using StatsBase, Plots

include("train_NN.jl")
include("utils_GTM.jl")
include("constants.jl")
result_dir = joinpath(root_path, "result/synthetic")
model_type = "GTM"
d = 5
intensity = 1
logranges = [  1.0 15.0
                0.1 10.0
                1.0 15.0
                0.01 10.0
                0.1  400.0
             ]

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

    param_names = [:x1, :x2, :x3, :x4, :x5]
    col_names = [:loss; param_names]
    
    # 初始化 DataFrame，指定列名和列类型（例如 Float64）
    history = DataFrame([Symbol(n) => Float64[] for n in col_names]...)

    best_loss = Inf
    wait = 0
    temp_loss = Inf
    if length(param)==0
        prior = Product(Uniform.(logranges[1:d, 1], logranges[1:d, 2]))
        param = rand(prior)
    end
    history = DataFrame()
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
        temp_loss_copy = exp.(-deepcopy(temp_loss))
        param_copy = deepcopy(param)

        row = NamedTuple{Tuple(col_names)}((temp_loss_copy, param_copy...))
        push!(history, row)
    end
    # history = calculate_burst(history)
    pred = pred_pdf_infe(model, param, 0:length(hist_yobs))
    temp_loss = hellinger2(pred, hist_yobs)
    return param, temp_loss,history
end
             
function loss_hellinger_map_baye(x::AbstractVector, model,hist_yobs, tt=tt)
    bufs = zeros(Threads.nthreads())
    Threads.@threads for i in 1:length(tt)
        ps = x
        # pred = pred_pdf(model, ps, 0:length(hist_yobs))
        pred = pred_pdf_infe(model, ps, 0:length(hist_yobs))
        bufs[Threads.threadid()] += hellinger2(pred, hist_yobs)
    end
    global AllSolutions
    global AllLoss
    push!(AllSolutions, deepcopy(x))
    loss_result =sum(bufs)
    push!(AllLoss, loss_result)
    return loss_result
end

function calculate_burst(data)
    data[:, :tauOn] = data[!,1]./data[!,2]
    data[:, :tauOff] = data[!,3]./data[!,4]
    data[:, :tauOnRatio] = data[:, :tauOn] ./(data[:, :tauOn]+data[:, :tauOff])
    data[:, :bf] = 1 ./ (data[:, :tauOn] .+ data[:, :tauOff])
    data[:,:bs] = data[:, :tauOn] .* data[!,5]
    return data
end

# 自定义带权 KDE 计算
function weighted_kde(x::Vector{Float64}, w::Vector{Float64}; bandwidth::Float64=0.3, gridsize::Int=200)
    xmin, xmax = minimum(x), maximum(x)
    xgrid = range(xmin, xmax, length=gridsize)
    density = zeros(Float64, length(xgrid))

    for (i, xi) in enumerate(xgrid)
        for (xj, wj) in zip(x, w)
            density[i] += wj * pdf(Normal(xj, bandwidth), xi)
        end
    end

    return xgrid, density
end

@load joinpath(DATA_DIR, "test_$model_type.jld2") params_arr_test sim_results_pro_test
@load joinpath(MODELWEIGHT_DIR, "model_stats_prob.jld2") model

estimates = []
s = SobolSeq(logranges[:,1] ,logranges[:,2])
ps_train = [ Sobol.next!(s) for i in 1:100 ]
start_point_num = 3
params_num = 100
params_name = ["kon","ron","koff","roff","mu"]
length_results_all = [] 
estimates_all = []

paramsArr_matrix = reduce(vcat,transpose.(params_arr_test[1:params_num]))
paramsArr_matrix_df = DataFrame(paramsArr_matrix,:auto)
CSV.write( joinpath(result_dir, "posteriorDist/true_params.csv"),paramsArr_matrix_df)
for i in eachindex(sim_results_pro_test[1:params_num])
    estimates = []
    AllSolutions = Vector{Float64}[]
    AllLoss = Vector{Float64}()
    hist_yobs = sim_results_pro_test[i]
    tt = [0,0]
    yobs = []
    length_results_i = Int[]
    p = Progress(1)
    Threads.@threads for j in 1:start_point_num
        opt_result = bboptimize(p -> loss_hellinger_map_baye(p, model,hist_yobs,tt),ps_train[j]; SearchRange = [ tuple(logranges[i,:]...) for i in 1:d ], TraceMode=:silent)
        push!(estimates, best_candidate(opt_result))
        # push!(length_results_i, length_result)
        ProgressMeter.next!(p)
    end
    # push!(length_results_all, length_results_i)
    push!(estimates_all, estimates)

    AllLoss=exp.(-AllLoss)
    AllSolutions_matrix = reduce(vcat,transpose.(AllSolutions))
    paramsArr_matrix = reduce(vcat,transpose.(params_arr_test[1:params_num]))

    AllSolution_df = DataFrame(AllSolutions_matrix,:auto)
    AllSolution_df[:, :tauOn] = AllSolution_df[!,1]./AllSolution_df[!,2]
    AllSolution_df[:, :tauOff] = AllSolution_df[!,3]./AllSolution_df[!,4]
    AllSolution_df[:, :tauOnRatio] = AllSolution_df[:, :tauOn] ./(AllSolution_df[:, :tauOn]+AllSolution_df[:, :tauOff])
    AllSolution_df[:, :bf] = 1 ./ (AllSolution_df[:, :tauOn] .+ AllSolution_df[:, :tauOff])
    AllSolution_df[:,:bs] = AllSolution_df[:, :tauOn] .* AllSolution_df[!,5]
    AllSolution_df[:,:loss] = AllLoss
    CSV.write( joinpath(result_dir, "posteriorDist/AllSolution_$i.csv"),AllSolution_df)
end

for i in 1:length(estimates_all)
    AllSolution_df = CSV.read(joinpath(result_dir, "posteriorDist/AllSolution_$i.csv"), DataFrame)

    for j in 1:length(estimates_all[i])
        initial_value = deepcopy(estimates_all[i][j])
        param, temp_loss, solution_nn= TX_inferrer(sim_results_pro_test[i],model;param=initial_value)
        solution_nn = calculate_burst(solution_nn)
        
        solution_nn = calculate_burst(solution_nn)
        AllSolution_df = vcat(AllSolution_df, solution_nn)
    end
    CSV.write( joinpath(result_dir, "posteriorDist/AllSolution_nn_$i.csv"),AllSolution_df)
end

param_index = 7
AllSolution_df_nn = CSV.read(joinpath(result_dir, "posteriorDist/AllSolution_nn_$param_index.csv"), DataFrame)
# # 模拟数据
# df = DataFrame(value = randn(1000), weight = rand(1000))
# 提取数据
x = AllSolution_df_nn[:,4]
w = StatsBase.normalize(AllSolution_df_nn.loss, 1) 
# 计算加权 KDE
xgrid, ydensity = weighted_kde(x, w)
# 绘图
plot(xgrid, ydensity, xlabel="Value", ylabel="Density", title="Weighted KDE (Manual)", legend=false)
vline!([params_arr_test[param_index][4]], color=:red, linewidth=2, linestyle=:dash, label="True Value")
savefig("D:\\academic_relate_code_two\\Nessie-main\\DeepTXcopy4\\DeepTX-main\\trainInferSSA\\result\\synthetic\\posteriorDist\\figure\\weighted_kde.png")




# 设置 2x2 子图布局
for param_index in 1:params_num

AllSolution_df_nn = CSV.read(joinpath(result_dir, "posteriorDist/AllSolution_nn_$param_index.csv"), DataFrame)

true_params = params_arr_test[param_index]
plt_layout = @layout [a b; c d]
plots_list = []


# 对每列进行加权 KDE 并画图
for i in 1:4
    
    x = AllSolution_df_nn[:, i]
    xgrid, ydensity = weighted_kde(x, w)
    
    p = plot(
        xgrid, ydensity,
        xlabel="Value", ylabel="Density",
        title="Param $i", legend=false
    )
    vline!(p, [true_params[i]], color=:red, linewidth=2, linestyle=:dash, label="True Value")
    push!(plots_list, p)
end

# 拼图
big_plot = plot(plots_list..., layout=plt_layout)
savefig(big_plot, "D:\\academic_relate_code_two\\Nessie-main\\DeepTXcopy4\\DeepTX-main\\trainInferSSA\\result\\synthetic\\posteriorDist\\figure\\weighted_kde_all$param_index.png")
end







history = DataFrame()
push!(history, (temp_loss, param...))

row = hcat(temp_loss, param')  # param' 是行向量
push!(history, vec(row))

param_names = [:param_1, :param_2, :param_3, :param_4, :param_5]
col_names = [:loss; param_names]
row = NamedTuple{Tuple(col_names)}((temp_loss, param...))
push!(history, row)

history = calculate_burst(history)

min_val = minimum(history[:, 3])
max_val = maximum(history[:, 3])

df_merged = vcat(AllSolution_df, history)

AllSolution_df = CSV.read(joinpath(result_dir, "posteriorDist/AllSolution_1.csv"), DataFrame)
AllSolution_df_sub = AllSolution_df[1:10128,:]
best_row = AllSolution_df_sub[argmin(AllSolution_df_sub.loss), 1:5]
values = collect(best_row)