using CSV
using DataFrames
using JLD2
include("../train_NN.jl")
include("contants.jl")
include("utils.jl")

root_path = "D:\\academic_relate_code_two\\Nessie-main\\momentsNN"
sim_type = "GTM"


@load joinpath(DATA_DIR, "train_$sim_type.jld2") params_arr_train sim_results_pro_train theo_stats_train

@load joinpath(DATA_DIR, "valid_$sim_type.jld2") params_arr_valid sim_results_pro_valid theo_stats_valid

params_arr_train = convert(Vector{Vector{Float32}},params_arr_train)
sim_results_pro_train = convert(Vector{Vector{Float32}},sim_results_pro_train)
theo_stats_train = convert(Vector{Vector{Float32}},theo_stats_train)

params_arr_valid = convert(Vector{Vector{Float32}},params_arr_valid)
sim_results_pro_valid = convert(Vector{Vector{Float32}},sim_results_pro_valid)
theo_stats_valid = convert(Vector{Vector{Float32}},theo_stats_valid)


data_step = 3000
data_end = 12001
params_arr_train_subs =  [params_arr_train[i:i+data_step] for i in collect(1:data_step:data_end)]
sim_results_pro_train_subs = [sim_results_pro_train[i:i+data_step] for i in collect(1:data_step:data_end)]
theo_stats_train_subs = [theo_stats_train[i:i+data_step] for i in collect(1:data_step:data_end)]


build_model(n_comps::Int, x::Int) = build_model(n_comps, [x])

""" Create a neural network with the specified number of components
    and hidden units per layer (we normally use 1 hidden layer) """
function build_model(n_comps::Int, x::Vector{Int}=[32])
    hidden_layers = [Dense(x[i-1], x[i], relu) for i in 2:length(x)]
    model = Chain(InputLayer(),
                  Dense(5, x[1], relu),
                  hidden_layers...,
                  MNBOutputLayer(x[end], n_comps)
            )
    MNBModel(model)
end


nsamples = 5
nunits = 2 .^(3:8)

models = []
@time for nunit in nunits
    
    _models = []
    
    for i in 1:nsamples
        # _train_data = (params_arr_train_subs[i],sim_results_pro_train_subs[i],theo_stats_train_subs[i])
        _train_data = (params_arr_train,sim_results_pro_train,theo_stats_train)

        _valid_data = (params_arr_valid[1:500],sim_results_pro_valid[1:500],theo_stats_valid[1:500])
        model = build_model(4, nunit)
        _, _ = train_NN!(model, _train_data, _valid_data; max_rounds=30, lr=0.01, batchsize=64)
        push!(_models, model)
    end
    
    push!(models, _models)

end

@save joinpath(weight_path, "models_nunits.jld2") models
@load joinpath(weight_path, "models_nunits.jld2") models
valid_data = (params_arr_valid[1:500],sim_results_pro_valid[1:500])
losses_hell = [ [ mean_loss(valid_data..., model; loss=loss_hellinger) for model in models[i] ] for i in 1:length(models)]

x = string.(nunits)
y = mean.(losses_hell)
ystd = std.(losses_hell)
ind = findmax(y)[2]
ymax = (y[ind] + ystd[ind]) * 1.01+0.2

plt_hell = plot(x, y, yerror=ystd, xlabel="number of neurons", ylabel="Hellinger distance",
                markerstrokecolor = "#0088c3ff", leg=false, grid=false, lw=1.5, ylim=(0., ymax),
                markerstrokewidth=1.5, tick_direction=:out, c=colorant"#0088c3ff",
                size = (290, 130), guidefontsize=6, tickfontsize=6, thickness_scaling=1.0,
                left_margin=-1Plots.mm, right_margin=-1Plots.mm, top_margin=-1Plots.mm, bottom_margin=0Plots.mm)

savefig(joinpath(figure_path, "hell_vs_neurons.pdf"))



nlayers = [[128], [128, 128], [64, 16], [128, 64, 32], [64, 32, 16], [64, 64, 32, 32], [16, 16, 16, 16]]
nsamples = 5

models = []
@time for nlayer in nlayers
    
    _models = []
    
    for i in 1:nsamples
        _train_data = (params_arr_train_subs[i],sim_results_pro_train_subs[i],theo_stats_train_subs[i])
        # _train_data = (params_arr_train,sim_results_pro_train,theo_stats_train)

        _valid_data = (params_arr_valid[1:500],sim_results_pro_valid[1:500],theo_stats_valid[1:500])
        model = build_model(4, nlayer)
        _, _ = train_NN!(model, _train_data, _valid_data; max_rounds=30, lr=0.01, batchsize=64)

        push!(_models, model)
    end
    
    push!(models, _models)

end

@save joinpath(weight_path, "models_nlayers.jld2") models
@load joinpath(weight_path, "models_nlayers.jld2") models

losses_hell = [ [ mean_loss(valid_data..., model; loss=loss_hellinger) for model in models[i] ] for i in 1:length(models)]
labels = ["128", "128-128", "64-16", "128-64-32", "64-32-16", "64-64-32-32", "16-16-16-16"]

plt_hell = plot(labels, mean.(losses_hell), yerror=std.(losses_hell), 
                xlabel="hidden layer architecture", ylabel="Hellinger distance", 
                markerstrokecolor = colorant"#0088c3ff", leg=false, grid=false, lw=1.5,
                markerstrokewidth=1.5, tick_direction=:out, c=colorant"#0088c3ff",
                size = (290, 130), guidefontsize=6, tickfontsize=6, thickness_scaling=1.0, xrotation = 30,
                left_margin=-1Plots.mm, right_margin=-2Plots.mm, top_margin=-1.5Plots.mm, bottom_margin=3Plots.mm)

savefig(joinpath(figure_path, "hell_vs_architecture.pdf"))

# ---------------------------------------------------------------------------------------------------
# Size of the dataset
# ---------------------------------------------------------------------------------------------------

nsamples = 5
nparams = [500,  1000, 5000, 10000,15000]

models = []
@time for nparam in nparams
    
    _models = []

    for i in 1:nsamples

        _train_data = (params_arr_train[1:nparams[i]],sim_results_pro_train[1:nparams[i]],theo_stats_train[1:nparams[i]])

        _valid_data = (params_arr_valid[1:500],sim_results_pro_valid[1:500],theo_stats_valid[1:500])
        
        model = build_model(4, 128)
        _, _ = train_NN!(model, _train_data, _valid_data; max_rounds=30, lr=0.01, batchsize=64)
        push!(_models, model)
    end
    
    push!(models, _models)

end

@save joinpath(weight_path, "models_nparams.jld2") models
@load joinpath(weight_path, "models_nparams.jld2") models

losses_hell = [ [ mean_loss(valid_data..., model; loss=loss_hellinger) for model in models[i] ] for i in 1:length(models)]

x = nparams
y = mean.(losses_hell)
ystd = std.(losses_hell)
ind = findmax(y)[2]
ymax = (y[ind] + ystd[ind]) * 1.01
xmax = x[end] * 1.01

plt_hell = plot(x, y, yerror=ystd, xlabel="dataset size", ylabel="Hellinger distance", 
                markerstrokecolor = colorant"#0088c3ff", leg=false, grid=false, lw=1.5, ylim=(0., ymax), xlim = (0., xmax),
                markerstrokewidth=1.5, tick_direction=:out, c=colorant"#0088c3ff",
                size = (290, 130), guidefontsize=6, tickfontsize=6, thickness_scaling=1.0,
                left_margin=-1Plots.mm, right_margin=0Plots.mm, top_margin=-1Plots.mm, bottom_margin=0Plots.mm)

savefig(joinpath(figure_path, "hell_vs_dataset.pdf"))


# ---------------------------------------------------------------------------------------------------
# Number of mixture components
# ---------------------------------------------------------------------------------------------------

nsamples = 5
ncomps = 1:10

models = []
@time for ncomp in ncomps
    
    _models = []
    
    for i in 1:nsamples
        _train_data = (params_arr_train_subs[i],sim_results_pro_train_subs[i],theo_stats_train_subs[i])

        _valid_data = (params_arr_valid[1:500],sim_results_pro_valid[1:500],theo_stats_valid[1:500])
        
        model = build_model(ncomp, 128)
        _, _ = train_NN!(model, _train_data, _valid_data; max_rounds=30, lr=0.01, batchsize=64)
        push!(_models, model)
    end
    
    push!(models, _models)

end

@save joinpath(weight_path, "models_ncomps.jld2") models
@load joinpath(weight_path, "models_ncomps.jld2") models

losses_hell = [ [ mean_loss(valid_data..., model; loss=loss_hellinger) for model in models[i] ] for i in 1:length(models)]

x = ncomps
y = mean.(losses_hell)
ystd = std.(losses_hell)
ind = findmax(y)[2]
ymax = (y[ind] + ystd[ind]) * 1.01

plt_hell = plot(x, y, yerror=ystd, xlabel="number of mixture components", ylabel="Hellinger distance", 
                markerstrokecolor = colorant"#0088c3ff", leg=false, grid=false, lw=1.5, ylim=(0., ymax),
                markerstrokewidth=1.5, tick_direction=:out, c=colorant"#0088c3ff", xticks=x,
                size = (290, 130), guidefontsize=6, tickfontsize=6, thickness_scaling=1.0,
                left_margin=-1Plots.mm, right_margin=-2Plots.mm, top_margin=-1Plots.mm, bottom_margin=0Plots.mm)

savefig(joinpath(figure_path, "hell_vs_ncomps.pdf"))
