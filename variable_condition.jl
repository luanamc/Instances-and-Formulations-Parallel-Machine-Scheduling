using Parameters

@with_kw mutable struct VariableCondition
    condition1::Array{Bool} = Array{Bool}(undef, 0)
    condition2::Array{Bool} = Array{Bool}(undef, 0)
end