using JuMP, Gurobi

function economic_dispatch(
    c::AbstractVector,
    a::AbstractMatrix,
    p_max::AbstractVector,
    p_min::AbstractVector,
    demand::Real,
)
    # Check that inputs are valid
    @assert length(c) == size(a, 2) == length(p_max) == length(p_min)
    @assert all(p_max .>= p_min)

    n = length(c)
    m = length(p_max)
    b = ones(m)

    # Define optimization problem using JuMP
    model = Model(with_optimizer(GLPK.Optimizer))
    @variable(model, p[1:m] >= 0)
    @objective(model, Min, dot(c, p))
    @constraint(model, a * p .== demand * b)
    @constraint(model, p .<= p_max)
    @constraint(model, p .>= p_min)
    optimize!(model)

    # Extract solution
    if termination_status(model) == MOI.OPTIMAL
        return value.(p)
    else
        error("No optimal solution found.")
    end
end