using JuMP, Gurobi

function economicdispatch(
    EDL::Matrix{Float64}, # Load
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
    GPini::Vector{Float64}, # generator initial output
    hydromap::Matrix{Int64}, # hydro map
    HAvail::Matrix{Float64}, # hydro availability
    renewablemap::Matrix{Int64}, # renewable map
    RAvail::Matrix{Float64}, # renewable availability
    U::Matrix{Int}; # Conventional generator status, 1 if on, 0 if off
    Horizon::Int = 1, # planning horizon
    Steps::Int = 12, # planning steps
    VOLL::Float64 = 9000.0, # value of lost load
    RM::Float64 = 0.2, # reserve margin
)::JuMP.Model
    ntimepoints = Horizon # number of time points
    nbus = size(EDL, 1) # number of buses
    ntrans = size(transmap, 1) # number of transmission lines
    nucgen = size(genmap, 1) # number of conventional generators
    nhydro = size(hydromap, 1) # number of hydro generators
    nrenewable = size(renewablemap, 1) # number of renewable generators

    # Define model
    edmodel = Model()

    # Define decision variables
    @variable(edmodel, f[1:ntrans, 1:ntimepoints]) # Transmission flow
    @variable(edmodel, theta[1:nbus, 1:ntimepoints]) # Voltage angle
    @variable(edmodel, guc[1:nucgen, 1:ntimepoints] >= 0) # Generator output
    @variable(edmodel, gh[1:nhydro, 1:ntimepoints] >= 0) # Hydro output
    @variable(edmodel, gr[1:nrenewable, 1:ntimepoints] >= 0) # Renewable output
    @variable(edmodel, s[1:nbus, 1:ntimepoints] >= 0) # Slack variable
    @variable(edmodel, θ[1:nbus, 1:ntimepoints]) # Phase angle 

    # Define objective function and constraints
    @objective(
        edmodel,
        Min,
        sum(GMC .* guc + GNLC .* U) / Steps +
        sum(VOLL .* s) / Steps
    )

    # Bus wise load balance constraints with transmission
    @constraint(
        edmodel,
        LoadBalance[z = 1:nbus, t = 1:ntimepoints],
        sum(genmap[:, z] .* guc[:, t]) +
        sum(hydromap[:, z] .* gh[:, t]) +
        sum(renewablemap[:, z] .* gr[:, t]) +
        sum(transmap[:, z] .* f[:, t]) +
        s[z, t] == EDL[z, t]
    )

    # Load balance constraints without transmission
    # @constraint(
    #     edmodel,
    #     LoadBalance[t = 1:ntimepoints],
    #     sum(guc[:, t]) + sum(gh[:, t]) + sum(gr[:, t]) + sum(s[:, t]) ==
    #     sum(EDL[:, t])
    # )

    # System reserve constraints
    # @constraint(
    #     edmodel,
    #     Reserve[t = 1:ntimepoints],
    #     sum(guc[:, t]) +
    #     sum(gh[:, t]) +
    #     sum(gr[:, t]) >= (1 + RM) * sum(EDL[:, t])
    # )

    # # Transmission capacity limits
    @constraint(
        edmodel,
        TXCapTo[l = 1:ntrans, t = 1:ntimepoints],
        f[l, t] <= TFmax[l]
    )
    @constraint(
        edmodel,
        TXCapFrom[l = 1:ntrans, t = 1:ntimepoints],
        f[l, t] >= -TFmax[l]
    )

    # # DCOPF constraints
    # set reference angle
    @constraint(edmodel, REFBUS[t = 1:ntimepoints], θ[1, t] == 0)
    @constraint(
        edmodel,
        DCOPTX[l = 1:ntrans, t = 1:ntimepoints],
        f[l, t] ==
        TX[l] * (
            θ[(findfirst(x -> x == 1, transmap[l, :])), t] -
            θ[(findfirst(x -> x == -1, transmap[l, :])), t]
        )
    )

    # Power flow constraints
    # @constraint(
    #     ucmodel,
    #     PowerFlow[i = 1:ntrans, t = 1:ntimep oints],
    #     f[i, t] ==
    #     (transmap[i, 1] - transmap[i, 2]) *
    #     TX[i] *
    #     (guc[transmap[i, 1], t] - guc[transmap[i, 2], t])
    # )

    # Conventional generator capacity limits
    @constraint(
        edmodel,
        UCCapU[i = 1:nucgen, t = 1:Horizon],
        guc[i, t] <= U[i] * GPmax[i]
    )
    @constraint(
        edmodel,
        UCCapL[i = 1:nucgen, t = 1:Horizon],
        guc[i, t] >= U[i] * GPmin[i]
    )
    # Hydro and renewable capacity limits
    @constraint(
        edmodel,
        HCap[i = 1:nhydro, t = 1:Horizon],
        gh[i, t] <= HAvail[i, t]
    )
    @constraint(
        edmodel,
        ReCap[i = 1:nrenewable, t = 1:Horizon],
        gr[i, t] <= RAvail[i, t]
    )
    # Ramping limits
    @constraint(
        edmodel,
        RUIni[i = 1:nucgen],
        guc[i, 1] - GPini[i] <= GRU[i] / Steps
    )
    @constraint(
        edmodel,
        RDIni[i = 1:nucgen],
        GPini[i] - guc[i, 1] <= GRD[i] / Steps
    )
    # if EDHorizon > 1
    # @constraint(
    #     ucmodel,
    #     RU[i = 1:nucgen, t = 1:Horizon-1],
    #     guc[i, t+1] - guc[i, t] <= GRU[i]
    # )
    # @constraint(
    #     ucmodel,
    #     RD[i = 1:nucgen, t = 1:Horizon-1],
    #     guc[i, t] - guc[i, t+1] <= GRD[i]
    # )
    return edmodel
end