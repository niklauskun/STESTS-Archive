using JuMP, DataFrames, CSV

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
    GPini::Vector{Float64},
    GMC::Vector{Float64},
    GNLC::Vector{Float64},
    GSUC::Vector{Float64},
    EDL::Matrix{Float64},
    EDHAvail::Matrix{Float64},
    EDRAvail::Matrix{Float64},
    ucmodel::JuMP.Model,
    edmodel::JuMP.Model;
    UCHorizon::Int = 24,
    EDHorizon::Int = 1,
    EDSteps::Int = 12,
    VOLL::Float64 = 9000.0,
)
    UCcost = Array{Float64}(undef, Nday)
    UCnetgen = Array{Float64}(undef, Nday, UCHorizon)
    UCgen = Array{Float64}(undef, Nday, UCHorizon)
    EDcost = Array{Float64}(undef, Nday * UCHorizon * EDSteps)

    set_optimizer_attribute(ucmodel, "MIPGap", 0.01)
    for d in 1:Nday
        # Update load and renewable generation
        LInput =
            convert(Matrix{Float64}, UCL[UCHorizon*(d-1)+1:UCHorizon*d, :]')
        HAvailInput =
            convert(Matrix{Float64}, HAvail[UCHorizon*(d-1)+1:UCHorizon*d, :]')
        RAvailInput =
            convert(Matrix{Float64}, RAvail[UCHorizon*(d-1)+1:UCHorizon*d, :]')
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
            UCnetgen[d, :] = sum(value.(ucmodel[:guc]), dims = 1)[:]
            UCgen[d, :] =
                sum(value.(ucmodel[:guc]), dims = 1)[:] +
                sum(value.(ucmodel[:gr]), dims = 1)[:] +
                sum(value.(ucmodel[:gh]), dims = 1)[:]
            U = value.(ucmodel[:u])
            V = value.(ucmodel[:v])
            W = value.(ucmodel[:w])
            S = value.(ucmodel[:s])
            Udf = DataFrame(U', :auto)
            Vdf = DataFrame(V', :auto)
            Wdf = DataFrame(W', :auto)
            Sdf = DataFrame(S', :auto)
            CSV.write("output/UCCommit.csv", Udf, append = true)
            CSV.write("output/UCStart.csv", Vdf, append = true)
            CSV.write("output/UCShut.csv", Wdf, append = true)
            CSV.write("output/UCSlack.csv", Sdf, append = true)
            UCGdf = DataFrame(value.(ucmodel[:guc])', :auto)
            CSV.write("output/UCGen.csv", UCGdf, append = true)
            # Update parameters in unit commitment model
            # GPini = value.(ucmodel[:guc])[:,UCHorizon]
            # Update initail generation output for next day TODO: do this after last ED
            # for i in axes(GPini, 1)
            #     set_normalized_rhs(ucmodel[:RUIni][i], GPini[i] + GRU[i])
            #     set_normalized_rhs(ucmodel[:RDIni][i], -GPini[i] + GRD[i])
            # end
            # must stay on and off constraints
            # GSU = fill(1, 71)
            # @constraint(
            #     ucmodel,
            #     MustOn[i = 1:71],
            #     sum(ucmodel[:u][i,1:GSU[i]]) == GSU[i]
            # )
            # Upscale unit commitment status for ED model
            EDV = zeros(size(V, 1), size(V, 2) * EDSteps)
            EDW = zeros(size(V, 1), size(V, 2) * EDSteps)
            for i in axes(V, 2)
                EDV[:, (i-1)*EDSteps+1] = V[:, i]
                EDW[:, (i-1)*EDSteps+1] = W[:, i]
            end

            for h in 1:UCHorizon
                for t in 1:EDSteps
                    ts = ((d - 1) * UCHorizon + h - 1) * EDSteps + t # time step
                    # Update status and strat up variable in objective function
                    @objective(
                        edmodel,
                        Min,
                        sum(
                            GMC .* edmodel[:guc] / EDSteps +
                            GNLC .* U[:, h] / EDSteps +
                            GSUC .* EDV[:, (h-1)*EDSteps+t],
                        ) +
                        sum(50 .* edmodel[:d] - 20 .* edmodel[:c]) / EDSteps +
                        sum(VOLL .* edmodel[:s]) / EDSteps
                    )
                    EDLInput =
                        convert(Matrix{Float64}, EDL[ts:ts+EDHorizon-1, :]')
                    EDHAvailInput = convert(
                        Matrix{Float64},
                        EDHAvail[ts:ts+EDHorizon-1, :]',
                    )
                    EDRAvailInput = convert(
                        Matrix{Float64},
                        EDRAvail[ts:ts+EDHorizon-1, :]',
                    )
                    for tp in 1:EDHorizon
                        # Ad hoc method to solve initial generation output for ED model
                        # if ts ==1 && tp == 1
                        if (h - 1) * EDSteps + t == 1
                            for i in axes(GPini, 1)
                                set_normalized_rhs(
                                    edmodel[:RUIni][i],
                                    GPini[i] +
                                    GRU[i] +
                                    (GPmin[i]) * EDV[i, (h-1)*EDSteps+t],
                                )
                                set_normalized_rhs(
                                    edmodel[:RDIni][i],
                                    -GPini[i] +
                                    GRD[i] +
                                    (GPmin[i]) * EDW[i, (h-1)*EDSteps+t],
                                )
                                # set_normalized_rhs(edmodel[:RUIni][i], GPini[i] + GRU[i])
                                # set_normalized_rhs(edmodel[:RDIni][i], -GPini[i] + GRD[i])
                            end
                        end
                        for z in axes(EDLInput, 1)
                            set_normalized_rhs(
                                edmodel[:LoadBalance][z, tp],
                                EDLInput[z, tp],
                            )
                        end
                        for i in axes(U, 1)
                            set_normalized_rhs(
                                edmodel[:UCCapU][i, tp],
                                U[i, h] * GPmax[i],
                            )
                            set_normalized_rhs(
                                edmodel[:UCCapL][i, tp],
                                U[i, h] * GPmin[i],
                            )
                        end
                        for i in axes(EDHAvailInput, 1)
                            set_normalized_rhs(
                                edmodel[:HCap][i, tp],
                                EDHAvailInput[i, tp],
                            )
                        end
                        for i in axes(EDRAvailInput, 1)
                            set_normalized_rhs(
                                edmodel[:ReCap][i, tp],
                                EDRAvailInput[i, tp],
                            )
                        end
                    end
                    # Solve economic dispatch model
                    optimize!(edmodel)
                    # # Extract solution
                    if termination_status(edmodel) == MOI.OPTIMAL
                        EDcost[ts] = objective_value(edmodel)
                        EDGPini = value.(edmodel[:guc])[:, 1]
                        EDGdf = DataFrame(EDGPini', :auto)
                        CSV.write("output/EDGen.csv", EDGdf, append = true)
                        EDSOCini = value.(edmodel[:e])[:, 1]
                        EDSOCinidf = DataFrame(EDSOCini', :auto)
                        CSV.write(
                            "output/EDSOCini.csv",
                            EDSOCinidf,
                            append = true,
                        )
                        EDS = value.(edmodel[:s])
                        EDSdf = DataFrame(EDS', :auto)
                        CSV.write("output/EDSlack.csv", EDSdf, append = true)
                        EDprice = dual.(edmodel[:LoadBalance])
                        EDpricedf = DataFrame(EDprice', :auto)
                        CSV.write(
                            "output/EDprice.csv",
                            EDpricedf,
                            append = true,
                        )
                        # Update initial generation output for next time step
                        for i in axes(GPini, 1)
                            if (h - 1) * EDSteps + t != UCHorizon * EDSteps
                                set_normalized_rhs(
                                    edmodel[:RUIni][i],
                                    EDGPini[i] +
                                    GRU[i] / EDSteps +
                                    (GPmin[i] + GRD[i] - GRU[i] / EDSteps) *
                                    EDV[i, (h-1)*EDSteps+t+1],
                                )
                                set_normalized_rhs(
                                    edmodel[:RDIni][i],
                                    -EDGPini[i] +
                                    GRD[i] / EDSteps +
                                    (EDGPini[i] - GRD[i] / EDSteps) *
                                    EDW[i, (h-1)*EDSteps+t+1],
                                )
                                # set_normalized_rhs(edmodel[:RUIni][i], EDGPini[i] + GRU[i] + GPmin[i]*EDV[i,(h - 1) * EDSteps + t + 1])
                                # set_normalized_rhs(edmodel[:RDIni][i], -EDGPini[i] + GRD[i] + GPmin[i]*EDW[i,(h - 1) * EDSteps + t + 1])
                            end
                        end
                        for i in axes(EDSOCini, 1)
                            set_normalized_rhs(
                                edmodel[:StorageSOCIni][i],
                                EDSOCini[i],
                            )
                        end
                    else
                        error(
                            "No optimal solution found for ED at hour $h step $t.",
                        )
                    end
                    # Update initial generation output for next day
                    if h == UCHorizon && t == EDSteps
                        for i in axes(GPini, 1)
                            set_normalized_rhs(
                                ucmodel[:RUIni][i],
                                EDGPini[i] + GRU[i],
                            )
                            set_normalized_rhs(
                                ucmodel[:RDIni][i],
                                -EDGPini[i] + GRD[i],
                            )
                        end
                    end
                end
            end

        else
            error("No optimal solution found for UC on day $d.")
        end
    end
    return UCcost, UCnetgen, UCgen, EDcost
end