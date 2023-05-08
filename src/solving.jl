using JuMP

# Solving unit commitment and economic dispatch model by iterativly updating parameters
function solving(Nday, D, ucmodel::JuMP.Model; Horizon::Int = 24)
    set_optimizer_attribute(ucmodel, "MIPGap", 0.0005)
    for d in 1:Nday
        # Update parameters
        DInput = convert(Matrix{Float64}, D[Horizon*(d-1)+1:Horizon*d, :]')
        nbus = size(DInput, 1)
        for t in 1:Horizon
            for b in 1:nbus
                set_normalized_rhs(ucmodel[:LoadBalance][b, t], DInput[b, t])
            end
        end
        # Solve unit commitment model
        optimize!(ucmodel)
        # Extract solution and solve economic dispatch model
        if termination_status(ucmodel) == MOI.OPTIMAL
            u = value.(ucmodel[:u])
            v = value.(ucmodel[:v])
            z = value.(ucmodel[:z])
            # Update parameters
            # ...
            # # Solve economic dispatch model
            optimize!(edmodel)
            # # Extract solution
            if termination_status(edmodel) == MOI.OPTIMAL
                # p = value.(p)
                # Update parameters
                # ...
            else
                error("No optimal solution found for ED.")
            end
        else
            error("No optimal solution found for UC.")
        end
    end
    return objective_value(ucmodel)
end