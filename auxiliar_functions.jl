include("instance.jl")
include("param.jl")

"""
Unroll instance for acess time indexed parameters: 
    machine set, job set, elegibility, precedence, processing_time, release date, deadline, period set
"""
function unroll_instance_time_indexed_formulation(instance::Instance)
    M = instance.M 
    J = instance.J
    # availability set
    M_j  = instance.M_j
    J_m  = instance.J_m
    prec_jj′ = instance.prec_jj′
    p_j  = instance.p_j
    r_j  = instance.r_j
    d_j  = instance.d_j
    T = instance.T
    Tmax_j = instance.Tmax_j

    return M, J, M_j, J_m, prec_jj′, p_j, r_j, d_j, T, Tmax_j
end

"""
Unroll instance for acess bucket indexed parameters: 
machine set, job set, elegibility, precedence, bucket set,k binary set, bucket size, processing_time, release date, deadline
"""
function unroll_instance_bucket_indexed_formulation(instance::Instance)
    M = instance.M 
    J = instance.J
    # availability set
    M_j  = instance.M_j
    J_m  = instance.J_m
    prec_jj′ = instance.prec_jj′
    B = instance.B
    K = instance.K
    Δ = instance.Δ
    P_j = instance.P_j
    π_j = instance.π_j
    R_j = instance.R_j
    ρ_j = instance.ρ_j
    D_j = instance.D_j
    δ_j = instance.δ_j
    
    return M, J, M_j, J_m, prec_jj′, B, K, Δ, P_j, π_j, R_j, ρ_j, D_j, δ_j 
end

"""
Unroll instance for acess instances parameters: 
all parameters
"""
function unroll_instance(instance::Instance)
    file_name = instance.file_name
    M = instance.M 
    J = instance.J
    # availability set
    M_j  = instance.M_j
    J_m  = instance.J_m
    prec_jj′ = instance.prec_jj′
    p_j  = instance.p_j
    r_j  = instance.r_j
    d_j  = instance.d_j
    L = instance.L
    L_m = instance.L_m
    T = instance.T
    Tmax_j = instance.Tmax_j
    B = instance.B
    K = instance.K
    Δ = instance.Δ
    P_j = instance.P_j
    π_j = instance.π_j
    R_j = instance.R_j
    ρ_j = instance.ρ_j
    D_j = instance.D_j
    δ_j = instance.δ_j

    return file_name, M, J, M_j, J_m, prec_jj′, p_j, r_j, d_j, L, L_m, T, Tmax_j, B, K, Δ, P_j, π_j, R_j, ρ_j, D_j, δ_j 
end

"""
Function to recalculate horizon if use release date and/or Deadline
"""
function recalculate_planning_horizon(param::Param)

    file_name, M, J, M_j, J_m, prec_jj′, p_j, r_j, d_j, L, L_m, T, Tmax_j, B, K, Δ, P_j, π_j, R_j, ρ_j, D_j, δ_j  = unroll_instance(param.instance) 

    if param.release_date && param.deadline
        start_time = minimum(r_j)
        end_time = maximum(d_j + p_j)
        param.instance.T = 1:(end_time - start_time + 1)
        param.instance.Tmax_j = [r_j[j]:d_j[j] for j in J]
        param.instance.r_j = r_j .- start_time .+ 1
        param.instance.d_j = d_j .- start_time .+ 1
        param.instance.R_j, param.instance.ρ_j = convert_parameter_to_bucket((r_j .- 1), Δ)
        param.instance.D_j, param.instance.δ_j = convert_parameter_to_bucket((d_j .- 1), Δ)
        param.instance.s_a = s_a .- start_time .+ 1
        param.instance.f_a = f_a .- start_time .+ 1
        param.instance.S_a, param.instance.θ_a = convert_parameter_to_bucket((s_a .- 1), Δ)
        param.instance.F_a, param.instance.ϕ_a = convert_parameter_to_bucket((f_a .- 1), Δ)
        param.instance.B = 1:maximum(D_j .+ P_j .+ 1)
    elseif param.release_date
        start_time = minimum(r_j)
        end_time = maximum(r_j) + sum(p_j) - 1
        param.instance.T = 1:(end_time - start_time + 1)
        param.instance.Tmax_j = [r_j[j]:((end_time - start_time + 1) - p_j[j] + 1) for j in J]
        param.instance.r_j = r_j .- start_time .+ 1
        param.instance.R_j, param.instance.ρ_j = convert_parameter_to_bucket((r_j .-1), Δ)
        param.instance.s_a = s_a .- start_time .+ 1
        param.instance.S_a, param.instance.θ_a = convert_parameter_to_bucket((s_a .- 1), Δ)
        param.instance.B = 1:maximum(R_j) .+ sum(P_j .+ 1)
    elseif param.deadline
        start_time = 1
        end_time = maximum(d_j .+ p_j .- 1)
        param.instance.T = 1:(end_time - start_time + 1)
        param.instance.Tmax_j = [1:d_j[j] for j in J]
        param.instance.B = 1:maximum(D_j .+ P_j .+ 1)
    end
end