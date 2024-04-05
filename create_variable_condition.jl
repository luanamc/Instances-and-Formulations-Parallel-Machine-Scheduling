include("params.jl")
include("auxiliar_functions.jl")


"""
Create Variable Condition
"""
function create_variable_condition(param::Param)

    model_name = param.model_name
    
    # Create a dictionary that maps model_name to the corresponding function
    model_functions = Dict(
        "TI" => create_variable_condition_time_indexed,
        "BI" => create_variable_condition_bucket_indexed
    )

    # Check if model_name exists in the dictionary, and if it does, call the associated function with keyword arguments
    if haskey(model_functions, model_name)
        return model_functions[model_name](param)
    else
        error("Unknown model_name: $model_name")
    end
end

"""
Time-Indexed Variable Condition
"""
function create_variable_condition_time_indexed(param::Param)

    # Set parameters
    instance, release_date, deadline = unroll_create_variable_condition(param)
    M, J, M_j, J_m, prec_jj′, p_j, r_j, d_j, T, Tmax_j = unroll_instance_time_indexed_formulation(instance)

    # Check conditions
    if release_date && deadline
        exists_x = [(t >= r_j[j] && t <= d_j[j] && t in Tmax_j[j] && j in J_m[m]) for m in M, j in J, t in T]
    elseif release_date
        exists_x = [(t >= r_j[j] && t in Tmax_j[j] && j in J_m[m]) for m in M, j in J, t in T]
    elseif deadline
        exists_x = [(t <= d_j[j] && t in Tmax_j[j] && j in J_m[m]) for m in M, j in J, t in T]
    else
        exists_x = [(t in Tmax_j[j] && j in J_m[m]) for m in M, j in J, t in T]
    end

    param.condition.condition1 = exists_x
    
end


"""
Bucket-Indexed Variable Condition
"""
function create_variable_condition_bucket_indexed(param::Param)

    # Set parameters
    instance, release_date, deadline = unroll_create_variable_condition(param)
    M, J, M_j, J_m, prec_jj′, B, K, Δ, P_j, π_j, R_j, ρ_j, D_j, δ_j = unroll_instance_bucket_indexed_formulation(instance)

    # Check conditions
    if release_date && deadline
        exists_u = [(b in B[max(1, R_j[j]):min(end - P_j[j] - k + 1, D_j[j])]) && (1.0 - π_j[j] < ρ_j[j] <= 1.0 || k==1 || b != R_j[j]) && (δ_j[j] <= 1.0 - π_j[j] || k==0 || b != D_j[j]) && (π_j[j] < 1.0 || k == 0) && (j in J_m[m]) for m in M, j in J, b in B, k in K]
        exists_v = [(b in B[max(P_j[j] + k, R_j[j] + P_j[j] + k - 1):min(end, D_j[j])]) && (1.0 - π_j[j] < ρ_j[j] <= 1.0 || k==1 || b != (R_j[j] + P_j[j] + k - 1)) && (δ_j[j] <= 1.0 - π_j[j] || k==0 || b != (D_j[j] + P_j[j] + k - 1)) && (π_j[j] < 1.0 || k == 0) && (j in J_m[m]) for m in M, j in J, b in B, k in K]
    elseif release_date
        exists_u = [(b in B[max(1, R_j[j]):end - P_j[j] - k + 1]) && (1.0 - π_j[j] < ρ_j[j] <= 1.0 || k==1 || b != R_j[j]) && (π_j[j] < 1.0 || k == 0) && (j in J_m[m]) for m in M, j in J, b in B, k in K]
        exists_v = [(b in B[max(P_j[j] + k, R_j[j] + P_j[j] + k - 1):end]) && (1.0 - π_j[j] < ρ_j[j] <= 1.0 || k==1 || b != (R_j[j] + P_j[j] + k - 1)) && (π_j[j] < 1.0 || k == 0) && (j in J_m[m]) for m in M, j in J, b in B, k in K]
    elseif deadline
        exists_u = [(b in B[1:min(end - P_j[j] - k + 1, D_j[j])]) && (δ_j[j] <= 1.0 - π_j[j] || k==0 || b != D_j[j]) && (π_j[j] < 1.0 || k == 0) && (j in J_m[m]) for m in M, j in J, b in B, k in K]
        exists_v = [(b in B[P_j[j] + k:min(end, D_j[j])]) && (δ_j[j] <= 1.0 - π_j[j] || k==0 || b != (D_j[j] + P_j[j] + k - 1)) && (π_j[j] < 1.0 || k == 0) && (j in J_m[m]) for m in M, j in J, b in B, k in K]
    else
        exists_u = [(b in B[1:end - P_j[j] - k + 1]) && (π_j[j] < 1.0 || k == 0) && (j in J_m[m]) for m in M, j in J, b in B, k in K]
        exists_v = [(b in B[P_j[j] + k:end]) && (π_j[j] < 1.0 || k == 0) && (j in J_m[m]) for m in M, j in J, b in B, k in K]
    end

    param.condition.condition1 = exists_u
    param.condition.condition2 = exists_v

end

"""
Unroll param for variable condition: 
instance, release date and deadline parameters
"""
function unroll_create_variable_condition(param::Param)
    
    instance = param.instance
    release_date = param.release_date
    deadline = param.deadline

    return instance, release_date, deadline
end