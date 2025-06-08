# Train GTM on training data 
using JLD2
using MAT
using Catalyst
using DataFrames
using CSV
include("./train_NN.jl")
include("constants.jl")
sim_type = "GTM"
const LOSSVALUE_DIR = joinpath(root_path, "result/lossVal")
max_rounds = 5

@load joinpath(DATA_DIR, "train_$sim_type.jld2") params_arr_train sim_results_pro_train theo_stats_train
params_arr_train = convert(Vector{Vector{Float32}},params_arr_train)
sim_results_pro_train = convert(Vector{Vector{Float32}},sim_results_pro_train)

@load joinpath(DATA_DIR, "valid_$sim_type.jld2") params_arr_valid sim_results_pro_valid theo_stats_valid
params_arr_valid = convert(Vector{Vector{Float32}},params_arr_valid)
sim_results_pro_valid = convert(Vector{Vector{Float32}},sim_results_pro_valid)

params_arr_valid = convert(Vector{Vector{Float32}},params_arr_valid)
sim_results_pro_valid = convert(Vector{Vector{Float32}},sim_results_pro_valid)

train_data = (params_arr_train,sim_results_pro_train,theo_stats_train)
valid_data = (params_arr_valid, sim_results_pro_valid, theo_stats_valid)


# Build the neural network
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

model = build_model(4, 128)

# Training

train_losses=Vector{Float32}()
valid_losses=Vector{Float32}()

# Three different loss functions, using moments only, probability density only, and moments plus probability density
# loss_type_li = ["only_stats", "only_prob","stats_prob"]
loss_type_li = ["stats_prob"]

for loss_type in loss_type_li
    if loss_type ==  "only_stats"
    model = build_model(4, 128)

    train_losses, valid_losses = train_NN!(model, train_data, valid_data;loss=loss_stats, max_rounds=max_rounds, lr=0.01, batchsize=64)
 
    @save joinpath(MODELWEIGHT_DIR, "model_$loss_type.jld2") model
    println(joinpath(MODELWEIGHT_DIR, "model_$loss_type.jld2"))
    train_losses_df = DataFrame([train_losses],:auto)
    CSV.write( joinpath(LOSSVALUE_DIR, "train_$loss_type.csv"),train_losses_df)

    elseif loss_type ==  "only_prob"
        model = build_model(4, 128)
        println("train_begin")
        println(joinpath(MODELWEIGHT_DIR, "model_$loss_type.jld2"))
        @time train_losses, valid_losses = train_NN!(model, train_data, valid_data;loss=loss_crossentropy, max_rounds=max_rounds, lr=0.01, batchsize=64)
        println("train_end")
        @save joinpath(MODELWEIGHT_DIR, "model_$loss_type.jld2") model
        println(joinpath(MODELWEIGHT_DIR, "model_$loss_type.jld2"))
        train_losses_df = DataFrame([train_losses],:auto)
        CSV.write( joinpath(LOSSVALUE_DIR, "train_$loss_type.csv"),train_losses_df)

    elseif loss_type ==  "stats_prob"
        model = build_model(4, 128)
        @time train_losses, valid_losses = train_NN!(model, train_data, valid_data;loss=loss_crossentropy_statical, max_rounds=max_rounds, lr=0.01, batchsize=64)
        @save joinpath(MODELWEIGHT_DIR, "model_$loss_type.jld2") model
        train_losses_df = DataFrame([train_losses],:auto)
        CSV.write( joinpath(LOSSVALUE_DIR, "train_$loss_type.csv"),train_losses_df)
    end
end

