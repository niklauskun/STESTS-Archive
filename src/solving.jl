using JuMP

# Solving unit commitment and economic dispatch model by iterativly updating parameters
function solving(Nday, D, ucmodel::JuMP.Model; Horizon::Int = 24)
    for d in 1:Nday
        # Update parameters
        DInput = convert(Matrix{Float64}, D[Horizon*(d-1)+1:Horizon*d, :]')
        for t in 1:Horizon
            set_normalized_rhs(ucmodel[:LoadBalance][t], sum(DInput[:, t]))
        end
        # Solve unit commitment model
        optimize!(ucmodel)
        # Extract solution
        if termination_status(ucmodel) == MOI.OPTIMAL
            # guc = value.(guc)
            # u = value.(u)
            # Update parameters
            # ...
            # # Solve economic dispatch model
            # optimize!(edmodel)
            # # Extract solution
            # if termination_status(edmodel) == MOI.OPTIMAL
            #     p = value.(p)
            #     # Update parameters
            #     # ...
            # else
            #     error("No optimal solution found.")
            # end
        else
            error("No optimal solution found.")
        end
    end
    return objective_value(ucmodel)
end