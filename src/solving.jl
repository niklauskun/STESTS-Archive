using JuMP

# Solving unit commitment and economic dispatch model by iterativly updating parameters
function solving(
    Nday::Int, 
    UCL::Matrix{Float64}, 
    HAvail::Matrix{Float64}, 
    RAvail::Matrix{Float64}, 
    GRU::Vector{Float64}, 
    GRD::Vector{Float64},
    GPmax::Vector{Float64},
    GPmin::Vector{Float64},
    EDL::Matrix{Float64},
    EDHAvail::Matrix{Float64},
    EDRAvail::Matrix{Float64},
    ucmodel::JuMP.Model,
    edmodel::JuMP.Model; 
    UCHorizon::Int = 24, 
    EDHorizon::Int = 1, 
    EDSteps::Int = 12
)
    UCcost = Array{Float64}(undef, Nday)
    UCnetgen = Array{Float64}(undef, Nday, UCHorizon)
    UCgen = Array{Float64}(undef, Nday, UCHorizon)
    EDcost = Array{Float64}(undef, Nday*UCHorizon*EDSteps)

    set_optimizer_attribute(ucmodel, "MIPGap", 0.01)
    for d in 1:Nday
        # Update load and renewable generation
        LInput = convert(Matrix{Float64}, UCL[UCHorizon*(d-1)+1:UCHorizon*d, :]')
        HAvailInput = convert(Matrix{Float64}, HAvail[UCHorizon*(d-1)+1:UCHorizon*d, :]')
        RAvailInput = convert(Matrix{Float64}, RAvail[UCHorizon*(d-1)+1:UCHorizon*d, :]')
        for h in 1:UCHorizon
            # set_normalized_rhs(ucmodel[:Reserve][h], (1 + 0.2) * sum(DInput[:, h])-sum(HAvailInput[:, h])-sum(RAvailInput[:, h]))
            for z in axes(LInput, 1)
                set_normalized_rhs(ucmodel[:LoadBalance][z, h], LInput[z, h])
            end
            for i in axes(HAvailInput, 1)
                set_normalized_rhs(ucmodel[:HCap][i, h], HAvailInput[i, h])
            end
            for i in axes(RAvailInput, 1)
                set_normalized_rhs(ucmodel[:ReCap][i, h], RAvailInput[i, h])
            end
        end
        # Solve unit commitment model
        optimize!(ucmodel)
        # Extract solution and solve economic dispatch model
        if termination_status(ucmodel) == MOI.OPTIMAL
            UCcost[d] = objective_value(ucmodel)
            UCnetgen[d,:] = sum(value.(ucmodel[:guc]), dims=1)[:]
            UCgen[d,:] = sum(value.(ucmodel[:guc]), dims=1)[:]+sum(value.(ucmodel[:gr]), dims=1)[:]+sum(value.(ucmodel[:gh]), dims=1)[:]
            U = value.(ucmodel[:u])
            # v = value.(ucmodel[:v])
            # w = value.(ucmodel[:w])
            # Update parameters in unit commitment model
            GPini = value.(ucmodel[:guc])[:,UCHorizon]
            # Update initail generation output for next day
            for i in axes(GPini, 1)
                set_normalized_rhs(ucmodel[:RUIni][i], GPini[i] + GRU[i])
                set_normalized_rhs(ucmodel[:RDIni][i], -GPini[i] + GRD[i])
            end

            for h in 1:UCHorizon
                for t in 1:EDSteps
                    ts = ((d - 1) * UCHorizon + h - 1) * EDSteps + t # time step
                    EDLInput = convert(Matrix{Float64}, EDL[ts:ts+EDHorizon-1, :]')
                    EDHAvailInput = convert(Matrix{Float64}, EDHAvail[ts:ts+EDHorizon-1, :]')
                    EDRAvailInput = convert(Matrix{Float64}, EDRAvail[ts:ts+EDHorizon-1, :]')
                    for tp in 1:EDHorizon
                        for z in axes(EDLInput, 1)
                            set_normalized_rhs(edmodel[:LoadBalance][z, tp], EDLInput[z, tp])
                        end
                        for i in axes(U, 1)
                            set_normalized_rhs(edmodel[:UCCapU][i, tp], U[i,h] * GPmax[i])
                            set_normalized_rhs(edmodel[:UCCapL][i, tp], U[i,h] * GPmin[i])
                        end
                        for i in axes(EDHAvailInput, 1)
                            set_normalized_rhs(edmodel[:HCap][i, tp], EDHAvailInput[i, tp])
                        end
                        for i in axes(EDRAvailInput, 1)
                            set_normalized_rhs(edmodel[:ReCap][i, tp], EDRAvailInput[i, tp])
                        end
                    end
                    # Solve economic dispatch model
                    optimize!(edmodel)
                    # # Extract solution
                    if termination_status(edmodel) == MOI.OPTIMAL
                        EDcost[ts] = objective_value(edmodel)
                        EDGPini = value.(edmodel[:guc])[:,1]
                        # Update initial generation output for next time step
                        for i in axes(GPini, 1)
                            set_normalized_rhs(edmodel[:RUIni][i], EDGPini[i] + GRU[i]/EDSteps)
                            set_normalized_rhs(edmodel[:RDIni][i], -EDGPini[i] + GRD[i]/EDSteps)
                        end
                    else
                        error("No optimal solution found for ED.")
                    end
                end
            end

        else
            error("No optimal solution found for UC.")
        end
    end
    return UCcost, UCnetgen, UCgen, EDcost
end