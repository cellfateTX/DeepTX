using JLD2
using Interp1d
using Flux
using Flux: crossentropy, onecold, onehotbatch, params, train!
using CSV
using DataFrames
using MAT
include("../train_NN.jl")
include("../utils.jl")
include("../GTMSSA.jl")
include("../constants.jl")
# root_path = "/public/home/zhjiajun/hzw/academic_code/DeepGTM/"

function simulatin_es_parameters(estimated_df,save_dir,seq_depth)
    if !isdir(save_dir)
        # 如果文件夹不存在，则创建文件夹
        mkdir(save_dir)
    end
    for i in 1:size(estimated_df)[1]
        # for i in 1:1000
            prob_result_array =  Vector{Float64}[]
            simPra = Vector(estimated_df[i,:])  
            push!(simPra,1.0)
            result = gillespieGTM(simPra,i)
            try
                insert_result = interp_data(result.data[:,3],result.time)
                proResult = calcuPro(insert_result)
                proResult = Hierarchical_model(proResult,seq_depth)
                push!(prob_result_array,proResult)
                save_prob_result(prob_result_array,joinpath(save_dir, "distribution_$i.csv"))
            catch e
                print(e)
                continue
            end
        end
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

#---------------verify params of scRNA-seq's distribution  inferred by model and its gillespie result

idu_path = "doseTen_estimated_model_stats_prob.csv"
# idu_save_path =joinpath(RESULT_DIR, "doseFiftySSADistribution/")

result_dir = joinpath(RESULT_DIR, "Fifty5FU/")
idu_data_path = joinpath(result_dir,"doseFifty_estimated_$MODEL_TYPE.csv")
idu_save_path =joinpath(result_dir, "doseFiftySSADistribution/")
estimated_df = DataFrame(CSV.File(joinpath(idu_data_path)))[:,1:5]
simulatin_es_parameters(estimated_df,idu_save_path,SEQ_DEPTH)

dmso_data_path = joinpath(result_dir,"doseOne_estimated_$MODEL_TYPE.csv")
dmso_save_path =joinpath(result_dir, "doseOneSSADistribution/")
estimated_df = DataFrame(CSV.File(joinpath(dmso_data_path)))[:,1:5]
simulatin_es_parameters(estimated_df,dmso_save_path,SEQ_DEPTH)

