using JuMP

# Solving unit commitment and economic dispatch model by iterativly updating parameters
function solving(Nday, D, ucmodel::JuMP.Model; Horizon::Int = 24)
    set_optimizer_attribute(ucmodel, "MIPGap", 0.0005)
    for d in 1:Nday
        # Update parameters
        DInput = convert(Matrix{Float64}, D[Horizon*(d-1)+1:Horizon*d, :]')
        nbus = size(DInput, 1)
        for t in 1:Horizon
            for z in 1:nbus
                set_normalized_rhs(ucmodel[:LoadBalance][z, t], DInput[z, t])
            end
        end
        # Solve unit commitment model
        optimize!(ucmodel)
        # Extract solution and solve economic dispatch model
        if termination_status(ucmodel) == MOI.OPTIMAL
            u = value.(ucmodel[:u])
            v = value.(ucmodel[:v])
            z = value.(ucmodel[:z])
            uccost = objective_value(ucmodel)
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
                error("No optimal solution found.")
            end
        else
            error("No optimal solution found.")
        end
    end
    return objective_value(ucmodel)
end