using JuMP
"""
    unitcommitment()::JuMP.Model

Read system parameters and construct unit commitment model.
Nik's note: In julia vectorized functions are not requireed for performance.
            For readability, I write constriants in for loops.

# Methods
unitcommitment(L::Matrix{Float64})::JuMP.Model -> Multi-bus ucmodel
unitcommitment(L::Vector{Float64})::JuMP.Model -> Single-bus ucmodel
"""
function unitcommitment(
    UCL::Matrix{Float64}, # Load
    genmap::Matrix{Int64}, # generator map
    GPmax::Vector{Float64}, # generator maximum output
    GPmin::Vector{Float64}, # generator minimum output
    GMC::Vector{Float64}, # generator marginal cost
    transmap::Matrix{Int64}, # transmission map
    TX::Vector{Float64}, # transmission reactance
    TFmax::Vector{Float64}, # transmission maximum flow,
    GNLC::Vector{Float64}, # generator non-load-carrying cost
    GRU::Vector{Float64}, # generator ramp up rate
    GRD::Vector{Float64}, # generator ramp down rate
    GSUC::Vector{Float64}, # generator start-up cost
    GUT::Vector{Int64}, # generator minimum up time
    GDT::Vector{Int64}, # generator minimum down time
    GPini::Vector{Float64}, # generator initial output
    UInput::Vector{Int64}, # Conventional generator status, 1 if on, 0 if off
    hydromap::Matrix{Int64}, # hydro map
    HAvail::Matrix{Float64}, # hydro availability
    renewablemap::Matrix{Int64}, # renewable map
    RAvail::Matrix{Float64}, # renewable availability
    storagemap::Matrix{Int64}, # storage map
    EPC::Vector{Float64}, # storage charging capacity
    EPD::Vector{Float64}, # storage discharging capacity
    Eeta::Vector{Float64}, # storage efficiency
    ESOC::Vector{Float64}, # storage state of charge capacity
    ESOCini::Vector{Float64}; # storage initial state of charge
    Horizon::Int = 24, # planning horizon
    VOLL::Float64 = 1000.0, # value of lost load
    RM::Float64 = 0.03, # reserve margin
)::JuMP.Model
    ntimepoints = Horizon # number of time points
    nbus = size(UCL, 1) # number of buses
    ntrans = size(transmap, 1) # number of transmission lines
    nucgen = size(genmap, 1) # number of conventional generators
    nhydro = size(hydromap, 1) # number of hydro generators
    nrenewable = size(renewablemap, 1) # number of renewable generators
    nstorage = size(storagemap, 1) # number of storage units

    # Define model
    ucmodel = Model()

    # Define decision variables
    @variable(ucmodel, f[1:ntrans, 1:ntimepoints]) # Transmission flow
    @variable(ucmodel, θ[1:nbus, 1:ntimepoints]) # Phase angle 
    @variable(ucmodel, guc[1:nucgen, 1:ntimepoints] >= 0) # Generator output
    @variable(ucmodel, u[1:nucgen, 1:ntimepoints], Bin) # Conventional generator status, 1 if on, 0 if off
    @variable(ucmodel, v[1:nucgen, 1:ntimepoints], Bin) # Conventional generator start-up decision, 1 if start-up, 0 otherwise
    @variable(ucmodel, w[1:nucgen, 1:ntimepoints], Bin) # Conventional generator shut-down decision, 1 if shut-down, 0 otherwise
    @variable(ucmodel, gh[1:nhydro, 1:ntimepoints] >= 0) # Hydro output
    @variable(ucmodel, gr[1:nrenewable, 1:ntimepoints] >= 0) # Renewable output
    @variable(ucmodel, s[1:nbus, 1:ntimepoints] >= 0) # Slack variable
    @variable(ucmodel, c[1:nstorage, 1:ntimepoints] >= 0) # Storage charging
    @variable(ucmodel, d[1:nstorage, 1:ntimepoints] >= 0) # Storage discharging
    @variable(ucmodel, e[1:nstorage, 1:ntimepoints] >= 0) # Storage energy level
    @variable(ucmodel, grr[1:nucgen, 1:ntimepoints] >= 0) # Conventional generator reserve

    # Define objective function and constraints
    @objective(
        ucmodel,
        Min,
        sum(GMC .* guc + GNLC .* u + GSUC .* v) +
        sum(50 .* d - 20 .* c) +
        sum(VOLL .* s)
    )

    # Bus wise load balance constraints with transmission
    @constraint(
        ucmodel,
        LoadBalance[z = 1:nbus, h = 1:ntimepoints],
        sum(genmap[:, z] .* guc[:, h]) +
        sum(hydromap[:, z] .* gh[:, h]) +
        sum(renewablemap[:, z] .* gr[:, h]) +
        sum(storagemap[:, z] .* d[:, h]) - sum(storagemap[:, z] .* c[:, h]) +
        sum(transmap[:, z] .* f[:, h]) +
        s[z, h] == UCL[z, h]
    )

    # Load balance constraints without transmission
    # @constraint(
    #     ucmodel,
    #     LoadBalance[h = 1:ntimepoints],
    #     sum(guc[:, h]) + sum(gh[:, h]) + sum(gr[:, h]) + sum(s[:, h]) ==
    #     sum(UCL[:, h])
    # )

    # System reserve constraints
    @constraint(
        ucmodel,
        UnitReserve1[i = 1:nucgen, h = 1:ntimepoints],
        grr[i, h] <= GPmax[i] * u[i, h] - guc[i, h]
    )

    @constraint(
        ucmodel,
        UnitReserve2[i = 1:nucgen, h = 1:ntimepoints],
        grr[i, h] <=  GRU[i]
    )

    @constraint(
        ucmodel,
        Reserve[h = 1:ntimepoints],
        sum(grr[:, h]) >= RM * maximum(sum(UCL, dims = 2))
    )

    # # Transmission capacity limits
    @constraint(
        ucmodel,
        TXCapTo[l = 1:ntrans, h = 1:ntimepoints],
        f[l, h] <= TFmax[l]
    )
    @constraint(
        ucmodel,
        TXCapFrom[l = 1:ntrans, h = 1:ntimepoints],
        f[l, h] >= -TFmax[l]
    )

    # # DCOPF constraints
    # set reference angle
    @constraint(ucmodel, REFBUS[h = 1:ntimepoints], θ[1, h] == 0)
    @constraint(
        ucmodel,
        DCOPTX[l = 1:ntrans, h = 1:ntimepoints],
        f[l, h] ==
        TX[l] * (
            θ[(findfirst(x -> x == 1, transmap[l, :])), h] -
            θ[(findfirst(x -> x == -1, transmap[l, :])), h]
        )
    )

    # Storage charge and discharge constraints
    @constraint(
        ucmodel,
        StorageCharge[i = 1:nstorage, h = 1:ntimepoints],
        c[i, h] <= EPC[i]
    )

    @constraint(
        ucmodel,
        StorageDischarge[i = 1:nstorage, h = 1:ntimepoints],
        d[i, h] <= EPD[i]
    )

    # Storage energy level constraints
    @constraint(
        ucmodel,
        StorageSOCCap[i = 1:nstorage, h = 1:ntimepoints],
        e[i, h] <= ESOC[i]
    )

    # Storage SOC evolution constraints
    @constraint(
        ucmodel,
        StorageSOCIni[i = 1:nstorage],
        e[i, 1] == ESOCini[i] + c[i, 1] * Eeta[i] - d[i, 1] / Eeta[i]
    )

    @constraint(
        ucmodel,
        StorageSOC[i = 1:nstorage, h = 2:ntimepoints],
        e[i, h] == e[i, h-1] + c[i, h] * Eeta[i] - d[i, h] / Eeta[i]
    )

    # Conventional generator capacity limits
    @constraint(
        ucmodel,
        UCCapU[i = 1:nucgen, h = 1:Horizon],
        guc[i, h] <= u[i, h] * GPmax[i]
    )
    @constraint(
        ucmodel,
        UCCapL[i = 1:nucgen, h = 1:Horizon],
        guc[i, h] >= u[i, h] * GPmin[i]
    )
    # Hydro and renewable capacity limits
    @constraint(
        ucmodel,
        HCap[i = 1:nhydro, h = 1:Horizon],
        gh[i, h] <= HAvail[i, h]
    )
    @constraint(
        ucmodel,
        ReCap[i = 1:nrenewable, h = 1:Horizon],
        gr[i, h] <= RAvail[i, h]
    )
    # Ramping limits
    # @constraint(ucmodel, RUIni[i = 1:nucgen], guc[i, 1] - GPini[i] <= GRU[i])
    # @constraint(ucmodel, RDIni[i = 1:nucgen], GPini[i] - guc[i, 1] <= GRD[i])
    # @constraint(
    #     ucmodel,
    #     RU[i = 1:nucgen, h = 1:Horizon-1],
    #     guc[i, h+1] - guc[i, h] <= GRU[i]
    # )
    # @constraint(
    #     ucmodel,
    #     RD[i = 1:nucgen, h = 1:Horizon-1],
    #     guc[i, h] - guc[i, h+1] <= GRD[i]
    # )
    @constraint(
        ucmodel,
        RUIni[i = 1:nucgen],
        guc[i, 1] - GPini[i] <= GRU[i] + GPmin[i] * v[i, 1]
    )
    @constraint(
        ucmodel,
        RDIni[i = 1:nucgen],
        GPini[i] - guc[i, 1] <= GRD[i] + GPmin[i] * w[i, 1]
    )
    @constraint(
        ucmodel,
        RU[i = 1:nucgen, h = 1:Horizon-1],
        guc[i, h+1] - guc[i, h] <= GRU[i] + GPmin[i] * v[i, h+1]
    )
    @constraint(
        ucmodel,
        RD[i = 1:nucgen, h = 1:Horizon-1],
        guc[i, h] - guc[i, h+1] <= GRD[i] + GPmin[i] * w[i, h+1]
    )

    # State transition constraints
    @constraint(
        ucmodel,
        ST0[i = 1:nucgen],
        u[i, 1] - UInput[i] == v[i, 1] - w[i, 1]
    )

    @constraint(
        ucmodel,
        ST1[i = 1:nucgen, h = 2:Horizon],
        u[i, h] - u[i, h-1] == v[i, h] - w[i, h]
    )

    # @constraint(
    #     ucmodel,
    #     ST1[i = 1:nucgen, h = 1:Horizon-1],
    #     u[i, h+1] - u[i, h] == v[i, h] - w[i, h]
    # )

    # @constraint(
    #     ucmodel,
    #     ST2[i = 1:nucgen, h = 1:Horizon-1],
    #      v[i, h] + w[i, h] <= 1
    # )

    # Minimum up/down time constraints
    @constraint(
        ucmodel,
        UTime[i = 1:nucgen, h = 1:Horizon],
        sum(v[i, max(1, h - GUT[i] + 1):h]) <= u[i, h]
    )
    @constraint(
        ucmodel,
        DTime[i = 1:nucgen, h = 1:Horizon],
        sum(w[i, max(1, h - GDT[i] + 1):h]) <= 1 - u[i, h]
    )

    return ucmodel
end