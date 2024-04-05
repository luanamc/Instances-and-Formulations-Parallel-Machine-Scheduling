using Parameters
using JuMP

# Include self files
include("instance.jl")
include("variable_condition.jl")

@with_kw mutable struct Param
    instance::Instance = Instance()
    model_name::String = ""
    model_objective::String = "min-completion"
    time_limit::Float64 = 3600.0
    condition::VariableCondition = VariableCondition()
    model::Model = Model()
    construction_time::Float64  = 0.0
    solve_time::Float64  = 0.0
    release_date::Bool = false
    deadline::Bool = false
    precedence::Bool = false
    availability::Bool = false
    log_file_name::String = ""
end
