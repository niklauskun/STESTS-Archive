using JuMP, DataFrames, CSV, Gurobi, Statistics

function getLastOneIndex(
    matrix::Array{Int,2},
    vector::Array{Int,1},
    vector2::Array{Int,1},
)
    m, n = size(matrix)
    result = Vector{Int}(undef, m)

    for i in 1:m
        idx = findlast(==(1), matrix[i, :])
        result[i] =
            idx === nothing ? max(0, vector[m] - 24) :
            max(0, vector2[m] - (24 + 1 - idx))
    end

    return result
end

function enforce_strictly_decreasing_vector(v, delta = 1e-5)
    indices = eachindex(v)
    for i in first(indices)+1:last(indices)
        if v[i] >= v[i-1]
            v[i-1] = v[i] + delta
        end
    end
end

# Solving unit commitment and economic dispatch model by iterativly updating parameters
function solving(
    params::STESTS.ModelParams,
    Nday::Int,
    strategic::Bool,
    DADBids::Matrix{Int64},
    DACBids::Matrix{Int64},
    RTDBids::Matrix{Int64},
    RTCBids::Matrix{Int64},
    ucmodel::JuMP.Model,
    ucpmodel::JuMP.Model,
    edmodel::JuMP.Model,
    bidmodels::Vector{Any},
    output_folder::String,
    PriceCap::Array{Float64};
    ESSeg::Int = 1,
    UCHorizon::Int = 24,
    EDHorizon::Int = 1,
    EDSteps::Int = 12,
    VOLL::Float64 = 9000.0,
    RM::Float64 = 0.03,
    ErrorAdjustment::Float64 = 1.0,
    LoadAdjustment::Float64 = 1.0,
)
    GSMC = repeat(params.GSMC, outer = (1, 1, UCHorizon))
    EDGSMC = repeat(params.GSMC, outer = (1, 1, EDHorizon))
    SU = zeros(Int, size(params.GPIni, 1)) # initial generator must on time
    SD = zeros(Int, size(params.GPIni, 1)) # initial generator down on time
    storagezone = [findmax(row)[2] for row in eachrow(params.storagemap)]
    segment_length = 5 ÷ ESSeg

    UCcost = Array{Float64}(undef, Nday)
    GMCcost = Array{Float64}(undef, Nday)
    GSMCcost = Array{Float64}(undef, Nday)
    VOLLcost = Array{Float64}(undef, Nday)
    EScost = Array{Float64}(undef, Nday)
    EDcost = Array{Float64}(undef, Nday * 24 * EDSteps)
    EDGMCcost = Array{Float64}(undef, Nday * 24 * EDSteps)
    EDGSMCcost = Array{Float64}(undef, Nday * 24 * EDSteps)
    EDVOLL = Array{Float64}(undef, Nday * 24 * EDSteps)
    EDEScost = Array{Float64}(undef, Nday * 24 * EDSteps)
    EDGPIni = Array{Float64}(undef, size(params.GPIni, 1))
    GPIniInput = params.GPIni
    all_UCprices_df = DataFrame()
    all_EDprices_df = DataFrame()
    vda = Array{Float64}(
        undef,
        size(params.storagemap, 1),
        ESSeg,
        24 * EDSteps + 1,
    )
    db = Array{Float64}(undef, size(params.storagemap, 1), ESSeg)
    cb = Array{Float64}(undef, size(params.storagemap, 1), ESSeg)
    AdjustedUCL = params.UCL * LoadAdjustment
    UCL_repeated = repeat(params.UCL, inner = (12, 1))
    LoadError = UCL_repeated - params.EDL
    AdjustedEDL = (UCL_repeated - LoadError * ErrorAdjustment) * LoadAdjustment
    Solar_repeated = repeat(params.SAvail, inner = (12, 1))
    SolarError = Solar_repeated - params.EDSAvail
    AdjustedEDSolar = Solar_repeated - SolarError * ErrorAdjustment
    Wind_repeated = repeat(params.WAvail, inner = (12, 1))
    WindError = Wind_repeated - params.EDWAvail
    AdjustedEDWind = Wind_repeated - WindError * ErrorAdjustment

    set_optimizer_attribute(ucmodel, "MIPGap", 0.001)
    for d in 1:Nday
        # Update load and renewable generation
        LInput = convert(
            Matrix{Float64},
            AdjustedUCL[24*(d-1)+1:24*d+UCHorizon-24, :]',
        )
        HAvailInput = convert(
            Matrix{Float64},
            params.HAvail[24*(d-1)+1:24*d+UCHorizon-24, :]',
        )
        # RAvailInput = convert(
        #     Matrix{Float64},
        #     params.RAvail[24*(d-1)+1:24*d+UCHorizon, :]',
        # )
        SAvailInput = convert(
            Matrix{Float64},
            params.SAvail[24*(d-1)+1:24*d+UCHorizon-24, :]',
        )
        WAvailInput = convert(
            Matrix{Float64},
            params.WAvail[24*(d-1)+1:24*d+UCHorizon-24, :]',
        )
        DADBidsInput =
            convert(Matrix{Float64}, DADBids[:, 24*(d-1)+1:24*d+UCHorizon-24])
        DACBidsInput =
            convert(Matrix{Float64}, DACBids[:, 24*(d-1)+1:24*d+UCHorizon-24])
        for h in 1:UCHorizon
            # set_normalized_rhs(ucmodel[:Reserve][h], (1 + 0.2) * sum(DInput[:, h])-sum(HAvailInput[:, h])-sum(RAvailInput[:, h]))
            for z in axes(LInput, 1)
                set_normalized_rhs(ucmodel[:LoadBalance][z, h], LInput[z, h])
                set_normalized_rhs(
                    ucmodel[:Reserve][h],
                    RM * sum(LInput, dims = 1)[h],
                )
                set_normalized_rhs(ucpmodel[:LoadBalance][z, h], LInput[z, h])
                set_normalized_rhs(
                    ucpmodel[:Reserve][h],
                    RM * sum(LInput, dims = 1)[h],
                )
            end
            for i in axes(HAvailInput, 1)
                set_normalized_rhs(ucmodel[:HCap][i, h], HAvailInput[i, h])
                set_normalized_rhs(ucpmodel[:HCap][i, h], HAvailInput[i, h])
            end
            # for i in axes(RAvailInput, 1)
            #     set_normalized_rhs(ucmodel[:ReCap][i, h], RAvailInput[i, h])
            #     set_normalized_rhs(ucpmodel[:ReCap][i, h], RAvailInput[i, h])
            # end
            for z in axes(SAvailInput, 1)
                set_normalized_rhs(ucmodel[:SCap][z, h], SAvailInput[z, h])
                set_normalized_rhs(ucpmodel[:SCap][z, h], SAvailInput[z, h])
            end
            for z in axes(WAvailInput, 1)
                set_normalized_rhs(ucmodel[:WCap][z, h], WAvailInput[z, h])
                set_normalized_rhs(ucpmodel[:WCap][z, h], WAvailInput[z, h])
            end
            for i in axes(DADBids, 1)
                set_objective_coefficient(
                    ucmodel,
                    ucmodel[:d][i, h],
                    DADBidsInput[i, h],
                )
                set_objective_coefficient(
                    ucmodel,
                    ucmodel[:c][i, h],
                    -DACBidsInput[i, h],
                )
                set_objective_coefficient(
                    ucpmodel,
                    ucpmodel[:d][i, h],
                    DADBidsInput[i, h],
                )
                set_objective_coefficient(
                    ucpmodel,
                    ucpmodel[:c][i, h],
                    -DACBidsInput[i, h],
                )
            end
        end

        # Solve unit commitment model
        @info "Solving day $d UC."
        optimize!(ucmodel)

        # Extract solution and solve economic dispatch model
        if termination_status(ucmodel) == MOI.OPTIMAL
            GMCcost[d] = sum(value.(ucmodel[:guc])[:, 1:24] .* params.GMC)
            GSMCcost[d] =
                sum(value.(ucmodel[:gucs])[:, :, 1:24] .* GSMC[:, :, 1:24])
            VOLLcost[d] = sum(value.(ucmodel[:s])[1:24] .* VOLL)
            # todo!!! update bids
            EScost[d] = sum(
                value.(ucmodel[:d])[:, 1:24] .* DADBids[:, 1:24] -
                value.(ucmodel[:c])[:, 1:24] .* DACBids[:, 1:24],
            )
            UCcost[d] =
                GMCcost[d] +
                GSMCcost[d] +
                VOLLcost[d] +
                EScost[d] +
                sum(
                    value.(ucmodel[:u])[:, 1:24] .* params.GNLC +
                    value.(ucmodel[:v])[:, 1:24] .* params.GSUC,
                )
            # UCcost[d] = objective_value(ucmodel)
            U = convert(Matrix{Int64}, round.(value.(ucmodel[:u])))
            V = convert(Matrix{Int64}, round.(value.(ucmodel[:v])))
            W = convert(Matrix{Int64}, round.(value.(ucmodel[:w])))
            S = value.(ucmodel[:s])
            Udf = DataFrame(U', :auto)
            Vdf = DataFrame(V', :auto)
            Wdf = DataFrame(W', :auto)
            Sdf = DataFrame(S', :auto)
            CSV.write(
                joinpath(output_folder, "UCCommit.csv"),
                Udf,
                append = true,
            )
            CSV.write(
                joinpath(output_folder, "UCStart.csv"),
                Vdf,
                append = true,
            )
            CSV.write(joinpath(output_folder, "UCShut.csv"), Wdf, append = true)
            CSV.write(
                joinpath(output_folder, "UCSlack.csv"),
                Sdf,
                append = true,
            )
            UCGdf = DataFrame(value.(ucmodel[:guc])', :auto)
            CSV.write(
                joinpath(output_folder, "UCGen.csv"),
                UCGdf,
                append = true,
            )
            # UCRedf = DataFrame(value.(ucmodel[:gr])', :auto)
            # CSV.write(
            #     joinpath(output_folder, "UCRenewable.csv"),
            #     UCRedf,
            #     append = true,
            # )
            UCSolardf = DataFrame(value.(ucmodel[:gs])', :auto)
            CSV.write(
                joinpath(output_folder, "UCSolar.csv"),
                UCSolardf,
                append = true,
            )
            UCWinddf = DataFrame(value.(ucmodel[:gw])', :auto)
            CSV.write(
                joinpath(output_folder, "UCWind.csv"),
                UCWinddf,
                append = true,
            )
            UCHdf = DataFrame(value.(ucmodel[:gh])', :auto)
            CSV.write(
                joinpath(output_folder, "UCHydro.csv"),
                UCHdf,
                append = true,
            )
            UCGZone = DataFrame(value.(ucmodel[:guc])' * params.genmap, :auto)
            CSV.write(
                joinpath(output_folder, "UCGenZone.csv"),
                UCGZone,
                append = true,
            )
            UCHydroZone =
                DataFrame(value.(ucmodel[:gh])' * params.hydromap, :auto)
            CSV.write(
                joinpath(output_folder, "UCHydroZone.csv"),
                UCHydroZone,
                append = true,
            )
            # UCRenewableZone =
            #     DataFrame(value.(ucmodel[:gr])' * params.renewablemap, :auto)
            # CSV.write(
            #     joinpath(output_folder, "UCRenewableZone.csv"),
            #     UCRenewableZone,
            #     append = true,
            # )
            UCSoCZone =
                DataFrame(value.(ucmodel[:e])' * params.storagemap, :auto)
            CSV.write(
                joinpath(output_folder, "UCSoCZone.csv"),
                UCSoCZone,
                append = true,
            )
            UCFdf = DataFrame(value.(ucmodel[:f])', :auto)
            CSV.write(
                joinpath(output_folder, "UCTrans.csv"),
                UCFdf,
                append = true,
            )

            # Update status and must up/down constraints for next day
            SU = getLastOneIndex(V[:, 1:24], SU, params.GUT)
            SD = getLastOneIndex(W[:, 1:24], SD, params.GDT)
            SUInt = zeros(Int, size(params.GPIni, 1), UCHorizon)
            SDInt = zeros(Int, size(params.GPIni, 1), UCHorizon)
            for i in axes(params.GPIni, 1)
                num_ones_SU = min(SU[i], UCHorizon)
                num_ones_SD = min(SD[i], UCHorizon)
                SUInt[i, 1:num_ones_SU] .= Int64(1)
                SDInt[i, 1:num_ones_SD] .= Int64(1)
            end

            for i in axes(params.GPIni, 1)
                set_normalized_rhs(ucmodel[:ST0][i], U[i, 24])
                set_normalized_rhs(ucmodel[:UTimeIni][i], SU[i])
                set_normalized_rhs(ucmodel[:DTimeIni][i], SD[i])
                for h in 1:UCHorizon
                    set_normalized_coefficient(
                        ucmodel[:UTimeIni][i],
                        ucmodel[:u][i, h],
                        SUInt[i, h],
                    )
                    set_normalized_coefficient(
                        ucmodel[:DTimeIni][i],
                        ucmodel[:u][i, h],
                        SDInt[i, h],
                    )
                end
            end

            # Calculate day-ahead price
            for i in axes(params.GPIni, 1)
                for h in 1:UCHorizon
                    set_normalized_rhs(
                        ucpmodel[:UnitReserve1][i, h],
                        params.GPmax[i] * U[i, h],
                    )
                    set_normalized_rhs(
                        ucpmodel[:UCGenSeg1][i, h],
                        params.GPmin[i] * U[i, h],
                    )
                    set_normalized_rhs(
                        ucpmodel[:UCCapU][i, h],
                        params.GPmax[i] * U[i, h],
                    )
                    set_normalized_rhs(
                        ucpmodel[:UCCapL][i, h],
                        params.GPmin[i] * U[i, h],
                    )
                end
                for h in 1:(UCHorizon-1)
                    set_normalized_rhs(
                        ucpmodel[:RU][i, h],
                        params.GRU[i] + params.GPmin[i] * V[i, h+1],
                    )
                    set_normalized_rhs(
                        ucpmodel[:RD][i, h],
                        params.GRD[i] + params.GPmin[i] * W[i, h+1],
                    )
                end
                set_normalized_rhs(
                    ucpmodel[:RUIni][i],
                    params.GRU[i] + params.GPmin[i] * V[i, 1] + GPIniInput[i],
                )
                set_normalized_rhs(
                    ucpmodel[:RDIni][i],
                    params.GRD[i] + params.GPmin[i] * W[i, 1] - GPIniInput[i],
                )
            end

            optimize!(ucpmodel)
            UCprice = dual.(ucpmodel[:LoadBalance][:, 1:24])
            UCpricedf = DataFrame(UCprice', :auto)
            CSV.write(
                joinpath(output_folder, "UCprice.csv"),
                UCpricedf,
                append = true,
            )
            all_UCprices_df = vcat(all_UCprices_df, UCpricedf)
            UCprice288 = repeat(UCprice', inner = (EDSteps, 1))

            # generate non-strategic bids
            for i in axes(RTDBids, 1)
                vda[i, :, :] = generate_value_function(
                    1 / EDSteps,
                    params.EPD[i] / params.ESOC[i] / EDSteps,
                    params.Eeta[i],
                    ESSeg,
                    UCprice288 * params.storagemap[i, :],
                )
            end

            # Solving economic dispatch model
            EDU = repeat(U, inner = (1, EDSteps))
            EDV = zeros(size(V, 1), size(V, 2) * EDSteps + EDHorizon)
            EDW = zeros(size(V, 1), size(V, 2) * EDSteps + EDHorizon)
            for i in axes(V, 2)
                EDU[:, (i-1)*EDSteps+1] = U[:, i]
                EDV[:, (i-1)*EDSteps+1] = V[:, i]
                EDW[:, (i-1)*EDSteps+1] = W[:, i]
            end

            for h in 1:24
                for t in 1:EDSteps
                    ts = ((d - 1) * 24 + h - 1) * EDSteps + t # time step
                    EDLInput = convert(
                        Matrix{Float64},
                        AdjustedEDL[ts:ts+EDHorizon-1, :]',
                    )
                    EDHAvailInput = convert(
                        Matrix{Float64},
                        params.EDHAvail[ts:ts+EDHorizon-1, :]',
                    )
                    # EDRAvailInput = convert(
                    #     Matrix{Float64},
                    #     params.EDRAvail[ts:ts+EDHorizon-1, :]',
                    # )
                    EDSAvailInput = convert(
                        Matrix{Float64},
                        AdjustedEDSolar[ts:ts+EDHorizon-1, :]',
                    )
                    EDWAvailInput = convert(
                        Matrix{Float64},
                        AdjustedEDWind[ts:ts+EDHorizon-1, :]',
                    )
                    EDDBidInput =
                        convert(Matrix{Float64}, RTDBids[:, ts:ts+EDHorizon-1])
                    EDCBidInput =
                        convert(Matrix{Float64}, RTCBids[:, ts:ts+EDHorizon-1])
                    if d > 1
                        last_24_UC = all_UCprices_df[(end-48+h):(end-25+h), :]
                        last_36_ED = all_EDprices_df[(end-35):end, :] .* EDSteps
                        predictors = vcat(last_24_UC, last_36_ED)
                    end
                    for tp in 1:EDHorizon
                        # Ad hoc method to solve initial generation output for ED model
                        if ts == 1 && tp == 1
                            # if h == 1 && t == 1 && tp == 1
                            for i in axes(params.GPIni, 1)
                                set_normalized_rhs(
                                    edmodel[:RUIni][i],
                                    params.GPIni[i] +
                                    params.GRU[i] +
                                    (params.GPmin[i]) *
                                    EDV[i, (h-1)*EDSteps+t+tp-1],
                                )
                                set_normalized_rhs(
                                    edmodel[:RDIni][i],
                                    -params.GPIni[i] +
                                    params.GRD[i] +
                                    (params.GPmin[i]) *
                                    EDW[i, (h-1)*EDSteps+t+tp-1],
                                )
                                # set_normalized_rhs(edmodel[:RUIni][i], GPini[i] + GRU[i])
                                # set_normalized_rhs(edmodel[:RDIni][i], -GPini[i] + GRD[i])
                            end
                        elseif h == 1 && t == 1 && tp == 1
                            for i in axes(params.GPIni, 1)
                                set_normalized_rhs(
                                    edmodel[:RUIni][i],
                                    EDGPIni[i] +
                                    params.GRU[i] / EDSteps +
                                    (params.GPmin[i]) * EDV[i, 1],
                                )
                                set_normalized_rhs(
                                    edmodel[:RDIni][i],
                                    -EDGPIni[i] +
                                    params.GRD[i] / EDSteps +
                                    (params.GPmax[i]) * EDW[i, 1],
                                )
                                # set_normalized_rhs(edmodel[:RUIni][i], GPini[i] + GRU[i])
                                # set_normalized_rhs(edmodel[:RDIni][i], -GPini[i] + GRD[i])
                            end
                            # TODO
                            # for i in axes(params.ESOCini, 1)
                            #     set_normalized_rhs(
                            #         edmodel[:StorageSOCIni][i],
                            #         params.ESOCini[i],
                            #     )
                            # end
                        end

                        # update energy storage bids
                        # update bids base on historical bids
                        for i in axes(RTDBids, 1)
                            db[i, :] .= EDDBidInput[i, tp] / EDSteps
                            cb[i, :] .= -EDCBidInput[i, tp] / EDSteps
                            # db[i, :] =
                            #     (
                            #         vda[i, :, (h-1)*EDSteps+t] ./
                            #         params.Eeta[i] .+ 10.0
                            #     ) / EDSteps
                            # cb[i, :] =
                            #     -vda[i, :, (h-1)*EDSteps+t] .* params.Eeta[i] /
                            #     EDSteps
                            # db[i, :] .=
                            #     (
                            #         params.storagemap[i, :]' * UCprice[:, h] /
                            #         0.9 + 10
                            #     ) / EDSteps
                            # cb[i, :] .=
                            #     -params.storagemap[i, :]' *
                            #     UCprice[:, h] *
                            #     0.9 / EDSteps
                            for s in 1:ESSeg
                                set_objective_coefficient(
                                    edmodel,
                                    edmodel[:d][i, s, tp],
                                    db[i, s],
                                )
                                set_objective_coefficient(
                                    edmodel,
                                    edmodel[:c][i, s, tp],
                                    cb[i, s],
                                )
                            end
                            cbdf = DataFrame(12 * cb[i, :]', :auto)
                            dbdf = DataFrame(12 * db[i, :]', :auto)
                            CSV.write(
                                joinpath(
                                    output_folder * "/NStrategic",
                                    "EDESCB*" * "$i" * ".csv",
                                ),
                                cbdf,
                                append = true,
                            )
                            CSV.write(
                                joinpath(
                                    output_folder * "/NStrategic",
                                    "EDESDB*" * "$i" * ".csv",
                                ),
                                dbdf,
                                append = true,
                            )
                        end
                        # update bids for strategic energy storage
                        if strategic == true
                            for (i, model) in bidmodels
                                if d == 1
                                    db[i, :] .= EDDBidInput[i, tp] / EDSteps
                                    cb[i, :] .= -EDCBidInput[i, tp] / EDSteps
                                    # db[i, :] =
                                    #     (
                                    #         vda[i, :, (h-1)*EDSteps+t] ./
                                    #         params.Eeta[i] .+ 10.0
                                    #     ) / EDSteps
                                    # cb[i, :] =
                                    #     -vda[i, :, (h-1)*EDSteps+t] .*
                                    #     params.Eeta[i] / EDSteps
                                    # db[i, :] .=
                                    #     (
                                    #         params.storagemap[i, :]' *
                                    #         UCprice[:, h] / 0.9 + 10
                                    #     ) / EDSteps
                                    # cb[i, :] .=
                                    #     -params.storagemap[i, :]' *
                                    #     UCprice[:, h] *
                                    #     0.9 / EDSteps
                                else
                                    v = model(
                                        predictors[
                                            :,
                                            findfirst(
                                                isequal(1),
                                                params.storagemap[i, :],
                                            ),
                                        ],
                                    )
                                    vavg = [
                                        mean(
                                            v[(i-1)*segment_length+1:i*segment_length],
                                        ) for i in 1:ESSeg
                                    ]
                                    enforce_strictly_decreasing_vector(vavg)
                                    cb[i, :] .= -vavg .* 0.9 ./ EDSteps
                                    db[i, :] .= (vavg ./ 0.9 .+ 10) ./ EDSteps
                                end
                                for s in 1:ESSeg
                                    set_objective_coefficient(
                                        edmodel,
                                        edmodel[:d][i, s, tp],
                                        db[i, s],
                                    )
                                    set_objective_coefficient(
                                        edmodel,
                                        edmodel[:c][i, s, tp],
                                        cb[i, s],
                                    )
                                end
                                # TODO: for single strategic storage now
                                cbdf = DataFrame(12 * cb[i, :]', :auto)
                                dbdf = DataFrame(12 * db[i, :]', :auto)
                                CSV.write(
                                    joinpath(
                                        output_folder * "/Strategic",
                                        "EDESCB*" * "$i" * ".csv",
                                    ),
                                    cbdf,
                                    append = true,
                                )
                                CSV.write(
                                    joinpath(
                                        output_folder * "/Strategic",
                                        "EDESDB*" * "$i" * ".csv",
                                    ),
                                    dbdf,
                                    append = true,
                                )
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
                                EDU[i, (h-1)*EDSteps+t+tp-1] * params.GPmax[i],
                            )
                            set_normalized_rhs(
                                edmodel[:UCCapL][i, tp],
                                EDU[i, (h-1)*EDSteps+t+tp-1] * params.GPmin[i],
                            )
                            set_normalized_rhs(
                                edmodel[:UCGenSeg1][i, tp],
                                EDU[i, (h-1)*EDSteps+t+tp-1] * params.GPmin[i],
                            )
                        end
                        for i in axes(EDHAvailInput, 1)
                            set_normalized_rhs(
                                edmodel[:HCap][i, tp],
                                EDHAvailInput[i, tp],
                            )
                        end
                        # for i in axes(EDRAvailInput, 1)
                        #     set_normalized_rhs(
                        #         edmodel[:ReCap][i, tp],
                        #         EDRAvailInput[i, tp],
                        #     )
                        # end
                        for z in axes(EDSAvailInput, 1)
                            set_normalized_rhs(
                                edmodel[:SCap][z, tp],
                                EDSAvailInput[z, tp],
                            )
                        end
                        for z in axes(EDWAvailInput, 1)
                            set_normalized_rhs(
                                edmodel[:WCap][z, tp],
                                EDWAvailInput[z, tp],
                            )
                        end
                        if EDHorizon > 1 && tp > 1
                            for i in axes(params.GPIni, 1)
                                set_normalized_rhs(
                                    edmodel[:RU][i, tp],
                                    params.GRU[i] / EDSteps +
                                    params.GPmin[i] *
                                    EDV[i, (h-1)*EDSteps+t+tp-1],
                                )
                            end
                            for i in axes(params.GPIni, 1)
                                set_normalized_rhs(
                                    edmodel[:RD][i, tp],
                                    params.GRD[i] / EDSteps +
                                    params.GPmax[i] *
                                    EDW[i, (h-1)*EDSteps+t+tp-1],
                                )
                            end
                        end
                    end
                    # Solve economic dispatch model
                    @info "Solving day $d hour $h step $t ED."
                    optimize!(edmodel)
                    # # Extract solution
                    if termination_status(edmodel) == MOI.OPTIMAL
                        EDGMCcost[ts] =
                            sum(value.(edmodel[:guc])[:, 1] .* params.GMC) /
                            EDSteps
                        EDGSMCcost[ts] =
                            sum(
                                value.(edmodel[:gucs])[:, :, 1] .*
                                EDGSMC[:, :, 1],
                            ) / EDSteps
                        EDVOLL[ts] =
                            sum(
                                value.(edmodel[:s])[:, :, 1] .*
                                PriceCap[:, :, 1],
                            ) / EDSteps
                        # todo!!! update bidmodel
                        # EDEScost[ts] =
                        #     sum(
                        #         value.(edmodel[:d])[:, :, 1] .*
                        #         RTDBids[:, (h-1)*EDSteps+t] -
                        #         value.(edmodel[:c])[:, :, 1] .*
                        #         RTCBids[:, (h-1)*EDSteps+t],
                        #     ) / EDSteps
                        EDEScost[ts] = sum(
                            value.(edmodel[:d])[:, :, 1] .*
                            (EDDBidInput[:, 1] / EDSteps) -
                            value.(edmodel[:c])[:, :, 1] .*
                            (EDCBidInput[:, 1] / EDSteps),
                        )
                        EDcost[ts] =
                            EDGMCcost[ts] +
                            EDGSMCcost[ts] +
                            EDVOLL[ts] +
                            EDEScost[ts] +
                            sum(
                                EDU[:, (h-1)*EDSteps+t] .* params.GNLC /
                                EDSteps +
                                EDV[:, (h-1)*EDSteps+t] .* params.GSUC,
                            )
                        EDGPIni = value.(edmodel[:guc])[:, 1]
                        EDGdf = DataFrame(EDGPIni', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDGen.csv"),
                            EDGdf,
                            append = true,
                        )
                        # EDReGen = value.(edmodel[:gr])[:, 1]
                        # EDReGendf = DataFrame(EDReGen', :auto)
                        # CSV.write(
                        #     joinpath(output_folder, "EDRenewable.csv"),
                        #     EDReGendf,
                        #     append = true,
                        # )
                        EDHydroGen = value.(edmodel[:gh])[:, 1]
                        EDHydroGendf = DataFrame(EDHydroGen', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDHydro.csv"),
                            EDHydroGendf,
                            append = true,
                        )
                        # CSV.write(
                        #     joinpath(output_folder, "EDDbid.csv"),
                        #     DataFrame(EDDBidInput[:, 1]', :auto),
                        #     append = true,
                        # )
                        # CSV.write(
                        #     joinpath(output_folder, "EDCbid.csv"),
                        #     DataFrame(EDCBidInput[:, 1]', :auto),
                        #     append = true,
                        # )
                        EDESD = value.(edmodel[:totald])[:, 1]
                        EDESDdf = DataFrame(EDESD', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDESD.csv"),
                            EDESDdf,
                            append = true,
                        )
                        EDESC = value.(edmodel[:totalc])[:, 1]
                        EDESCdf = DataFrame(EDESC', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDESC.csv"),
                            EDESCdf,
                            append = true,
                        )
                        EDSegSOCini = value.(edmodel[:e])[:, :, 1]
                        EDSOCini = value.(edmodel[:totale])[:, 1]
                        EDSOCinidf = DataFrame(EDSOCini', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDSOCini.csv"),
                            EDSOCinidf,
                            append = true,
                        )
                        EDTrans = value.(edmodel[:f])[:, 1]
                        EDTransdf = DataFrame(EDTrans', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDTrans.csv"),
                            EDTransdf,
                            append = true,
                        )
                        EDS = value.(edmodel[:totals])[:, 1]
                        EDSdf = DataFrame(EDS', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDSlack.csv"),
                            EDSdf,
                            append = true,
                        )
                        EDprice = dual.(edmodel[:LoadBalance])[:, 1]
                        EDpricedf = DataFrame(EDprice', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDprice.csv"),
                            EDpricedf,
                            append = true,
                        )
                        all_EDprices_df = vcat(all_EDprices_df, EDpricedf)
                        EDESRev =
                            EDprice[storagezone] .* (
                                value.(edmodel[:totald])[:, 1] -
                                value.(edmodel[:totalc])[:, 1]
                            )
                        CSV.write(
                            joinpath(output_folder, "EDESRev.csv"),
                            DataFrame(EDESRev', :auto),
                            append = true,
                        )
                        EDGZone = DataFrame(
                            value.(edmodel[:guc])[:, 1]' * params.genmap,
                            :auto,
                        )
                        CSV.write(
                            joinpath(output_folder, "EDGenZone.csv"),
                            EDGZone,
                            append = true,
                        )
                        EDHydroZone = DataFrame(
                            value.(edmodel[:gh])[:, 1]' * params.hydromap,
                            :auto,
                        )
                        CSV.write(
                            joinpath(output_folder, "EDHydroZone.csv"),
                            EDHydroZone,
                            append = true,
                        )
                        # EDRenewableZone = DataFrame(
                        #     value.(edmodel[:gr])[:, 1]' * params.renewablemap,
                        #     :auto,
                        # )
                        # CSV.write(
                        #     joinpath(output_folder, "EDRenewableZone.csv"),
                        #     EDRenewableZone,
                        #     append = true,
                        # )
                        EDSolarZone =
                            DataFrame(value.(edmodel[:gs])[:, 1]', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDSolarZone.csv"),
                            EDSolarZone,
                            append = true,
                        )
                        EDWindZone =
                            DataFrame(value.(edmodel[:gw])[:, 1]', :auto)
                        CSV.write(
                            joinpath(output_folder, "EDWindZone.csv"),
                            EDWindZone,
                            append = true,
                        )
                        EDESDisZone = DataFrame(
                            value.(edmodel[:totald])[:, 1]' * params.storagemap,
                            :auto,
                        )
                        CSV.write(
                            joinpath(output_folder, "EDESDisZone.csv"),
                            EDESDisZone,
                            append = true,
                        )
                        EDESChaZone = DataFrame(
                            value.(edmodel[:totalc])[:, 1]' * params.storagemap,
                            :auto,
                        )
                        CSV.write(
                            joinpath(output_folder, "EDESChaZone.csv"),
                            EDESChaZone,
                            append = true,
                        )
                        EDSOCZone = DataFrame(
                            value.(edmodel[:totale])[:, 1]' * params.storagemap,
                            :auto,
                        )
                        CSV.write(
                            joinpath(output_folder, "EDSoCZone.csv"),
                            EDSOCZone,
                            append = true,
                        )
                        # Update initial generation output for next time step
                        for i in axes(params.GPIni, 1)
                            # if (h - 1) * EDSteps + t != UCHorizon * EDSteps
                            set_normalized_rhs(
                                edmodel[:RUIni][i],
                                EDGPIni[i] +
                                params.GRU[i] / EDSteps +
                                (params.GPmin[i]) * EDV[i, (h-1)*EDSteps+t+1],
                            )
                            set_normalized_rhs(
                                edmodel[:RDIni][i],
                                -EDGPIni[i] +
                                params.GRD[i] / EDSteps +
                                (params.GPmax[i]) * EDW[i, (h-1)*EDSteps+t+1],
                            )
                            # set_normalized_rhs(edmodel[:RUIni][i], EDGPini[i] + GRU[i] + GPmin[i]*EDV[i,(h - 1) * EDSteps + t + 1])
                            # set_normalized_rhs(edmodel[:RDIni][i], -EDGPini[i] + GRD[i] + GPmin[i]*EDW[i,(h - 1) * EDSteps + t + 1])
                            # end
                        end
                        for i in axes(EDSOCini, 1)
                            for s in 1:ESSeg
                                set_normalized_rhs(
                                    edmodel[:StorageSOCIni][i, s],
                                    EDSegSOCini[i, s],
                                )
                            end
                        end
                    else
                        println("Model is infeasible. Starting IIS analysis...")
                        compute_conflict!(edmodel)
                        iis_model, _ = copy_conflict(edmodel)
                        open("diagnose_output_ed.txt", "w") do f
                            return write(f, string(iis_model))
                        end
                        print(iis_model)

                        error(
                            "No optimal solution found for ED at hour $h step $t.",
                        )
                    end
                    # Update initial generation output for next day
                    if h == 24 && t == 1
                        GPIniInput = EDGPIni
                        for i in axes(params.GPIni, 1)
                            set_normalized_rhs(
                                ucmodel[:RUIni][i],
                                EDGPIni[i] + params.GRU[i],
                            )
                            set_normalized_rhs(
                                ucmodel[:RDIni][i],
                                -EDGPIni[i] + params.GRD[i],
                            )
                            # set_normalized_rhs(
                            #     ucpmodel[:RUIni][i],
                            #     EDGPini[i] + GRU[i],
                            # )
                            # set_normalized_rhs(
                            #     ucpmodel[:RDIni][i],
                            #     -EDGPini[i] + GRD[i],
                            # )
                        end

                        # TODO: update initial storage SOC
                        for i in axes(EDSOCini, 1)
                            set_normalized_rhs(
                                ucmodel[:StorageSOCIni][i],
                                EDSOCini[i],
                            )
                            set_normalized_rhs(
                                ucpmodel[:StorageSOCIni][i],
                                EDSOCini[i],
                            )
                        end
                    end
                end
            end

        else
            println("Model is infeasible. Starting IIS analysis...")
            compute_conflict!(ucmodel)
            iis_model, _ = copy_conflict(ucmodel)
            open("diagnose_output.txt", "w") do f
                return write(f, string(iis_model))
            end
            error("No optimal solution found for UC on day $d.")
        end
    end
    EDcostdf = DataFrame(cost = EDcost)
    CSV.write(joinpath(output_folder, "EDcost.csv"), EDcostdf, append = true)
    return UCcost, EDcost
end