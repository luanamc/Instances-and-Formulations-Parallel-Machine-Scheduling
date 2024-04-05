# Get the directory of the current script or notebook
current_directory = @__DIR__ 
# Set the current working directory to the script's directory
cd(current_directory)

using Dates
using JuMP, Gurobi

# Include self files
include("params.jl")
include("auxiliar_functions.jl")
include("create_variable_condition.jl")
include("deterministic_models.jl")

println("------------ Start: "*param.instance.file_name*"  ------------")
println("Solution: "*param.model_name*" model")
println("Recalculating planning horizon")
if param.release_date || param.deadline
    recalculate_planning_horizon(param)
end

println("Creating variable condition")
condition = create_variable_condition(param)

println("Constructing model")
param.model, param.construction_time = generate_deterministic_model(param)
remain_time = max(param.time_limit - param.construction_time,0)

if remain_time > 0
    set_optimizer_attribute(param.model, "TimeLimit", remain_time)

    println("Optimizing model")
    param.solve_time = @elapsed begin
        optimize!(param.model)
else
    println("Optimization not called, because construction time exceedes the time limit.")
end
println("------------ End: "*param.instance.file_name*"  ------------")

