using Parameters

function fill_parameter_J_m(M::UnitRange{Int64}, J::UnitRange{Int64}, M_j::Vector{Vector{Int64}})
    J_m = [[] for _ in M]
    [append!(J_m[m],j) for j in J for m in M_j[j] if !(j in J_m[m])]
    return J_m
end

function convert_parameter_to_bucket(parameter::Array{Int64}, bucketSize::Int64)

    parameter_in_bucket = floor.(Int, parameter./bucketSize) .+ 1
    parameter_in_fraction_bucket = parameter_in_bucket .- parameter./bucketSize

    return parameter_in_bucket, parameter_in_fraction_bucket
end

# Create an  structure called Instance with all needed parameters to perform parallel machines scheduling problems with time-indexed and bucket-indexed formulation
@with_kw mutable struct Instance
    
    filename::String
    ## Problem Sets
    M::UnitRange{Int64}  # Set of machines
    J::UnitRange{Int64}  # Set of jobs
    A::UnitRange{Int64} # Set of availability
    ## Problem Parameters
    M_j::Vector{Vector{Int64}}  # Subset of machines can execute job j ∈ J
    p_j::Vector{Int64}  # Processing time of job j ∈ J
    r_j::Vector{Int64}  # Release date of job j ∈ J
    d_j::Vector{Int64}  # Deadline j ∈ J
    prec_jj′::Matrix{Int64} # 1, if job j precedes job j' and jobs j,j' ∈ J
    A_m::Vector{Vector{Int64}}  # Subset of availability in machine m ∈ M
    s_a::Vector{Int64} # Start time a
    f_a::Vector{Int64} # Finish time a
    ## Calculated parameters
    J_m::Vector{Vector{Int64}}  # Subset of jobs can execute m ∈ M
    ## Continuos Sets
    L::UnitRange{Int64}  # Set of position/location
    L_m::Vector{UnitRange{Int64}} 
    ## Time-Indexed Sets
    T::UnitRange{Int64}  # Set of periods
    ## Time-Indexed Parameters 
    Tmax_j::Vector{UnitRange{Int64}}  # Maximum number of periods per job j ∈ J
    #st_jj′          # setup time between jobs j,j' ∈ J
    ## Bucket-Indexed Sets
    B::UnitRange{Int64}  # Set of buckets
    K::UnitRange{Int64}  # Binary set
    ## Bucket-Indexed Parameters
    Δ::Int64 # Bucket size
    P_j::Vector{Int64}    # Processing time in bucket units of job j ∈ J, the minimum number of buckets occupied by job j 
    π_j::Vector{Float64}  # The proportion of a bucket necessary to indicate that the processing time of the job j does not occupy the entire bucket
    R_j::Vector{Int64}  # Release date in bucket units of job j ∈ J
    ρ_j::Vector{Float64}# The maximum proportion of a bucket which can be occupied in bucket Rj
    D_j::Vector{Int64}  # Deadline in bucket units of job j ∈ J
    δ_j::Vector{Float64}# The minimum proportion of a bucket which can be occupied in bucket Dj
    S_a::Vector{Int64}
    θ_a::Vector{Float64}
    F_a::Vector{Int64}
    ϕ_a::Vector{Float64}

    function Instance(filename::String, M::UnitRange{Int64}, J::UnitRange{Int64}, M_j::Vector{Vector{Int64}}, p_j::Vector{Int64}, r_j::Vector{Int64}, d_j::Vector{Int64}; prec_jj′::Matrix{Int64} = zeros(Int, 0, 0), A::UnitRange{Int64} = 1:0, A_m::Vector{Vector{Int64}} = [[]], s_a::Vector{Int64} = zeros(Int,0), f_a::Vector{Int64} = zeros(Int,0))
        J_m = fill_parameter_J_m(M, J, M_j)
        L = J 
        L_m = [L for _ in M]
        T = 1:sum(p_j)
        Tmax_j = [1:(sum(p_j) - p_j[j] + 1) for j in J]
        K = 0:1
        Δ = minimum(p_j)
        P_j, π_j = convert_parameter_to_bucket(p_j, Δ)
        R_j, ρ_j = convert_parameter_to_bucket((r_j .-1), Δ)
        D_j, δ_j = convert_parameter_to_bucket((d_j .-1), Δ)
        S_a, θ_a = convert_parameter_to_bucket((s_a .-1), Δ)
        F_a, ϕ_a = convert_parameter_to_bucket((f_a .-1), Δ)
        B = 1:sum(P_j .+ 1)
        new(filename, M, J, A, M_j, p_j, r_j, d_j, prec_jj′, A_m, s_a, f_a, J_m, L, L_m, T, Tmax_j, B, K, Δ, P_j, π_j, R_j, ρ_j, D_j, δ_j,S_a, θ_a,F_a, ϕ_a)
    end
end
