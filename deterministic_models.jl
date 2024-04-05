# Define packages
using JuMP, Gurobi
myGurobi = ()-> Gurobi.Optimizer(Gurobi.Env())

# Include self files
include("param.jl")
include("auxiliar_functions.jl")

"""
Generate deterministic model
"""
function generate_deterministic_model(param::Param)
    
    model_name = param.model_name
    
    # Create a dictionary that maps model_name to the corresponding function
    model_functions = Dict(
        "TI" => generate_time_index_formulation_deterministic_model,
        "BI" => generate_bucket_index_formulation_deterministic_model,
        "MTG" => generate_multiple_time_grid_formulation_deterministic_model
    )

    # Check if model_name exists in the dictionary, and if it does, call the associated function with keyword arguments
    if haskey(model_functions, model_name)
        return model_functions[model_name](param)
    else
        error("Unknown model_name: $model_name")
    end
end

"""
Multiple Time Grid Formulation
"""
function generate_multiple_time_grid_formulation_deterministic_model(param::Param)
    
    # Param
    instance, model_name, model_objective, condition, release_date, deadline = unroll_generate_deterministic_model(param)

    # Instance
    M, J, M_j, J_m, prec_jj′, p_j, r_j, d_j, T = unroll_instance_time_indexed_formulation(instance)

    # Initialize model
    construction_time = @elapsed begin
        
        model = Model(myGurobi)
        
        # Decision variables
        @variable(model, TS[M,T] >= 0)
        @variable(model, TF[M,T] >= 0)
        @variable(model, P[M,T] >= 0)
        # 1, if job j ∈ J the kth job in machine i ∈ I where k ∈ K; 0, otherwise
        @variable(model, x[M, J, T], Bin)
        @variable(model, x_free[M, T], Bin)

        if model_objective == "min-completion"
            # Minimize completion time
            @objective(model, Min, sum(TF[m,t] for m in M, t in T))
        elseif model_objective == "min-makespan"
            # value of makespan
            @variable(model, makespan >= 0)

            # Minimize completion time
            @objective(model, Min, makespan)

            # Constraint to calculate makespan
            @constraint(model, calculate_makespan[m in M], 
                makespan >= TF[m,maximum(T)] 
            )
        end

        # Constraints
        @constraint(model, nooverlap[m in M, t in T[1:end-1]],
            TS[m,t+1] >= TF[m,t]
        )

        @constraint(model, calculate_p[m in M, t in T],
            TF[m,t] - TS[m,t] == P[m,t]
        )

        @constraint(model, noorderbefore[m in M, t in T],
            sum(x[m,j,t] for j in J) + x_free[m,t] == 1
        )

        @constraint(model, noorderbefore_x_free[m in M, t in T[1:end-1]],
            x_free[m,t] <= x_free[m,t+1]
        )
        
        @constraint(model, onetime[j in J],
            sum(x[m,j,t] for m in M, t in T) == 1
        )

        @constraint(model, durationslottimegrid[m in M, t in T],
            P[m,t] == sum(p_j[j]*x[m,j,t] for j in J)
        )

        if release_date
            @constraint(model, release[m in M, t in T],
                TS[m,t] >= sum(r_j[j]*x[m,j,t] for j in J)
            )
        end

        if deadline
            @constraint(model, deadline[m in M, t in T],
                TF[m,t] <= sum((d_j[j]+p_j[j])*x[m,j,t] for j in J) + maximum(T)*x_free[m,t]
            )
        end

    end

    return model, construction_time  # Return the generated model

end

"""
Time-Indexed Formulation
"""
function generate_time_index_formulation_deterministic_model(param::Param)
    
    # Param
    instance, model_name, model_objective, condition, release_date, deadline = unroll_generate_deterministic_model(param)

    # Instance
    M, J, M_j, J_m, prec_jj′, p_j, r_j, d_j, T, Tmax_j = unroll_instance_time_indexed_formulation(instance)
    
    # get variable condition 
    exists_x = condition.condition1

    # Initialize model
    construction_time = @elapsed begin
        model = Model(myGurobi)

        # Decision variables
        # 1, if job j ∈ J starts in machine i ∈ I at a period t ∈ T; 0, otherwise
        @variable(model, x[m in M, j in J, t in T; exists_x[m,j,t]], Bin) 

        if model_objective == "min-completion"
            # Minimize completion time
            @objective(model, Min, sum((t + p_j[j] - 1)*x[m,j,t] for m in M, j in J, t in T if exists_x[m,j,t]))

        elseif model_objective == "min-makespan"
            # value of makespan
            @variable(model, makespan >= 0)

            # Minimize completion time
            @objective(model, Min, makespan)

            # Constraint to calculate makespan
            @constraint(model, calculate_makespan[j in J], 
                makespan >= sum((t + p_j[j] - 1)*x[m,j,t] for m in M, t in T if exists_x[m,j,t])
            )
        end

        # Constraints
        # Ensure that a job is executed at once
        @constraint(model, jobexecutedonce[j in J],
            sum(x[m,j,t] for m in M, t in T if exists_x[m,j,t]) == 1
        )
        
        # Ensure that a machine cannot perform more than one job per period
        @constraint(model, onejobperperiod[m in M, t in T],
            sum(x[m,j,t1] for j in J, t1 in max(1, t - p_j[j] + 1):t if exists_x[m,j,t1]) <= 1
        )
    end

    return model, construction_time  # Return the generated model

end

"""
Bucket-Indexed Formulation
"""
function generate_bucket_index_formulation_deterministic_model(param::Param)
    
    # Param
    instance, model_name, model_objective, condition, release_date, deadline = unroll_generate_deterministic_model(param)

    # Instance
    M, J, M_j, J_m, prec_jj′, B, K, Δ, P_j, π_j, R_j, ρ_j, D_j, δ_j = unroll_instance_bucket_indexed_formulation(instance)
    
    # get variable condition 
    exists_u = condition.condition1
    exists_v = condition.condition2

    # Initialize model
    construction_time = @elapsed begin
        
        model = Model(myGurobi)

        # Decision variables
        # 1, indicates if the machine m executes the job j at bucket b; 0, otherwise
        @variable(model, x[m in M, j in J, b in B, k in K; exists_u[m,j,b,k+1]], Bin)
        # indicates the proportion of bucket b is consumed by the job j
        @variable(model, 0.0 <= u[m in M, j in J, b in B, k in K; exists_u[m,j,b,k+1]] <= 1.0)
        # indicates the proportion of bucket b is consumed by the job j
        @expression(model, v[m in M, j in J, b in B, k in K; exists_v[m,j,b,k+1]], (2 - k - π_j[j])*x[m,j,b - P_j[j] - k + 1,k] - u[m,j,b - P_j[j] - k + 1,k])

        if model_objective == "min-completion"
            # Minimize completion time
            @objective(model, Min, sum((b + P_j[j] - π_j[j])*x[m,j,b,k] - u[m,j,b,k] for m in M, j in J, b in B, k in K if exists_u[m,j,b,k+1]))
        elseif model_objective == "min-makespan"
            # value of makespan
            @variable(model, makespan >= 0)

            # Minimize completion time
            @objective(model, Min, makespan)

            # Constraint to calculate makespan
            @constraint(model, calculate_makespan[j in J], 
                makespan >= sum((b + P_j[j] - π_j[j])*x[m,j,b,k] - u[m,j,b,k] for m in M, b in B, k in K if exists_u[m,j,b,k+1])
            )
        end
        

        # Constraints
        # Ensure that a job is executed at once
        @constraint(model, jobexecutedonce[j in J],
            sum(x[m,j,b,k] for m in M, b in B, k in K if exists_u[m,j,b,k+1]) == 1
        )

        # Ensure that a machine cannot start more than one job per bucket
        @constraint(model, onestartjobperbucket[m in M, b in B],
            sum(x[m,j,b1,k] for j in J, k in K, b1 = max(1, b - P_j[j] - k + 2):b if exists_u[m,j,b1,k+1]) <= 1
        )

        # Ensure that a machine cannot overlap jobs in a bucket
        @constraint(model, overlapjobs[m in M, b in B],
            sum((exists_u[m,j,b,k+1] ? u[m,j,b,k] : 0.0) + (exists_v[m,j,b,k+1] ? v[m,j,b,k] : 0.0) 
            + sum(x[m,j,b1,k] for b1 = max(1, b - P_j[j] - k + 2):(b - 1) if exists_u[m,j,b1,k+1]) for j in J, k in K) <= 1
        )
        
        # Lower bound variable u 
        @constraint(model, LBvaru[m in M, j in J, b in B, k in K; exists_u[m,j,b,k+1]],
            ((1 - k)*(1 - π_j[j]) + 1/Δ)*x[m,j,b,k] <= u[m,j,b,k]
        )

        # Upper bound variable u 
        @constraint(model, UBvaru[m in M, j in J, b in B, k in K; exists_u[m,j,b,k+1]],
            u[m,j,b,k] <= (1 - k*π_j[j])*x[m,j,b,k] 
        )

        if release_date
            # condition to generate constranints
            condition_release = [((b == R_j[j]) && (((k == 0) && (1.0 - π_j[j] < ρ_j[j] <= 1.0)) || ((k == 1) && (ρ_j[j] <= 1.0 - π_j[j])))) for m in M, j in J,  b in B, k in K]

            # Upper bound u variable considering release date
            @constraint(model, UBvarurelease[m in M, j in J, b = R_j[j], k in K; condition_release[m,j,b,k+1] && exists_u[m,j,b,k+1]],
                u[m,j,b,k] <= ρ_j[j]*x[m,j,b,k] 
            )
        end
        if deadline
            # condition to generate constranints
            condition_deadline = [((b == D_j[j]) && (((k == 0) && (1.0 - π_j[j] < δ_j[j] <= 1.0)) || (k == 1) && (δ_j[j] <= 1.0 - π_j[j]))) for m in M, j in J, b in B, k in K]

            # Lower bound u variable considering release date
            @constraint(model, UBvarudeadline[m in M, j in J, b = D_j[j], k in K; condition_deadline[m,j,b,k+1] && exists_u[m,j,b,k+1]],
                δ_j[j]*x[m,j,b,k] <= u[m,j,b,k]
            )
        end
    end

    return model, construction_time  # Return the generated model

end

"""
Unroll param for genereta deterministc model: 
instance, model name, model objective, condition, release date and deadline parameters
"""
function unroll_generate_deterministic_model(param::Param)
    
    instance = param.instance
    model_name = param.model_name
    model_objective = param.model_objective
    condition = param.condition
    release_date = param.release_date
    deadline= param.deadline

    return instance, model_name, model_objective, condition, release_date, deadline
end

