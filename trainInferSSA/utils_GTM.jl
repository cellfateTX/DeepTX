using MAT
using Interp1d

function saveMatFile(data,file_path,file_name)
    sim_param_results_df = DataFrame(data,:auto)
    file = matopen(file_path, "w")
    write(file, file_name, Matrix(sim_param_results_df))
    close(file)
end

function loadMatFile(file_path,var_name)
    
    file = matopen(file_path)
    var_names = keys(file)
    stat_the = read(file, var_name) # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    stat_the
end

function interp_data(yi::Vector{Int64},xi::Vector{Float64})
    yi = convert(Vector{Float64},yi)  
    x = collect(20000:0.05:30000)
    a=Vector{Float32}()
    for mode in INTERP_MODE_LIST
        if string(mode)=="Previous"
            b_index = 1
            f = interp(xi[b_index:length(xi)], yi[b_index:length(xi)], mode) # get an interpolation function
            y = f.(x) # Do interpolation
            a = y
        end
    end
    a
end

function calcuPro(data_counts)
    probArrNum = maximum(data_counts)
    probArrNum = convert(Int64,probArrNum)
    proArr = zeros(probArrNum+1)
    data_counts = convert(Vector{Int64},data_counts)
    for i in eachindex(data_counts)
        # num_index = convert(Int64,data_counts[i]+1)
       
        if data_counts[i]>=0
            num_index = data_counts[i]+1
            proArr[num_index]=proArr[num_index]+1
        end
    end
    proArr=proArr./length(data_counts)
    proArr
end

# transform the counts data array to probArr
function convertCountsToProb(data_counts)
    probArrNum = maximum(data_counts)
    probArrNum = convert(Int64, probArrNum)
    proArr = zeros(probArrNum + 1)
    data_counts = convert(Vector{Int64}, data_counts)
    for i in eachindex(data_counts)
        if data_counts[i] >= 0
            num_index = data_counts[i] + 1
            proArr[num_index] = proArr[num_index] + 1
        end
    end
    proArr = proArr ./ length(data_counts)
    proArr
end

# function linear_regression(X::AbstractVector,Y::AbstractVector;log_form=false)
#     if log_form
#         X=log.(X)
#         Y=log.(Y)
#     end
#     data = DataFrame(X=X, Y=Y)
# line = lm(@formula(Y ~ X), data)
# round(r2(line),digits=2) 
# end

function save_prob_result(data,file_path)
    df = DataFrame(data,:auto)
    CSV.write(file_path,df)
end

function convertCountToProb(data_counts)
    probArrNum = maximum(data_counts)
    probArrNum = convert(Int64,probArrNum)
    proArr = zeros(probArrNum+1)
    data_counts = convert(Vector{Int64},data_counts)
    for i in eachindex(data_counts)
        # num_index = convert(Int64,data_counts[i]+1)
        num_index = data_counts[i]+1
        proArr[num_index]=proArr[num_index]+1
    end
    proArr=proArr./length(data_counts)
    proArr
end

function slicematrix(A::AbstractMatrix)
    return [A[i, :] for i in 1:size(A,1)]
end

function hellinger2(aa, bb)
    ret = 0.0
    for (a, b) in zip(aa, bb)
        ret += (sqrt(a) - sqrt(b))^2
    end
    ret
end
