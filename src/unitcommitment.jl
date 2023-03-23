using JuMP
"""
    unitcommitment()::JuMP.Model

Read system parameters and construct unit commitment model.
Nik's note: In julia vectorized functions are not requireed for performance.
            For readability, I write constriants in for loops.

# Methods
unitcommitment(D::Matrix{Float64})::JuMP.Model -> Multi-bus ucmodel
unitcommitment(D::Vector{Float64})::JuMP.Model -> Single-bus ucmodel
"""
function unitcommitment(
    D::Matrix{Float64}, # demand
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
    GUT::Vector{Float64}, # generator minimum up time
    GDT::Vector{Float64}, # generator minimum down time
    GPIni::Vector{Float64}, # generator initial output
    hydromap::Matrix{Int64}, # hydro map
    HAvail::Matrix{Float64}, # hydro availability
    renewablemap::Matrix{Int64}, # renewable map
    RAvail::Matrix{Float64}; # renewable availability
    Horizon::Int = 24, # planning horizon
    VOLL::Float64 = 9000.0, # value of lost load
    RM::Float64 = 0.0, # reserve margin
)::JuMP.Model
    ntimepoints = Horizon # number of time points
    nbus = size(D, 2) # number of buses
    ntrans = size(transmap, 1) # number of transmission lines
    nucgen = size(genmap, 1) # number of conventional generators
    nhydro = size(hydromap, 1) # number of hydro generators
    nrenewable = size(renewablemap, 1) # number of renewable generators

    # Define model
    ucmodel = Model()

    # Define decision variables
    @variable(ucmodel, f[1:ntrans, 1:ntimepoints]) # Transmission flow
    @variable(ucmodel, theta[1:nbus, 1:ntimepoints]) # Voltage angle
    @variable(ucmodel, guc[1:nucgen, 1:ntimepoints] >= 0) # Generator output
    @variable(ucmodel, u[1:nucgen, 1:ntimepoints], Bin) # Conventional generator status, 1 if on, 0 if off
    @variable(ucmodel, v[1:nucgen, 1:ntimepoints], Bin) # Conventional generator start-up descison, 1 if start-up, 0 otherwise
    @variable(ucmodel, z[1:nucgen, 1:ntimepoints], Bin) # Conventional generator shut-down decision, 1 if shut-down, 0 otherwise
    @variable(ucmodel, gh[1:nhydro, 1:ntimepoints] >= 0) # Hydro output
    @variable(ucmodel, gr[1:nrenewable, 1:ntimepoints] >= 0) # Renewable output
    @variable(ucmodel, s[1:nbus, 1:ntimepoints] >= 0) # Slack variable

    # Define objective function and constraints
    @objective(
        ucmodel,
        Min,
        sum(GMC .* guc + GNLC .* u + GSUC .* v) + sum(VOLL .* s)
    )

    # Bus wise load balance constraints with transmission
    # @constraint(
    #     ucmodel,
    #     LoadBalance[i = 1:nbus, t = 1:ntimepoints],
    #     sum(genmap[:, i] .* guc[:, t]) +
    #     sum(hydromap[:, i] .* gh[:, t]) +
    #     sum(renewablemap[:, i] .* gr[:, t]) +
    #     sum(transmap[i, :] .== 1) .* f[transmap[i, :].==1, t] +
    #     s[i, t] == D[i, t]
    # )
    @constraint(
        ucmodel,
        LoadBalance[t = 1:ntimepoints],
        sum(guc[:, t]) + sum(gh[:, t]) + sum(gr[:, t]) + sum(s[:, t]) ==
        sum(D[:, t])
    )

    # set_name(LoadBalance, "LoadBalance")

    # System reserve constraints, individual bus reserve in current version
    # @constraint(
    #     ucmodel,
    #     Reserve[i = 1:nbus, t = 1:ntimepoints],
    #     sum(genmap[:, i] .* guc[:, t]) +
    #     sum(hydromap[:, i] .* gh[:, t]) +
    #     sum(renewablemap[:, i] .* gr[:, t]) >= (1 + RM) * D[t, i]
    # )

    # # Transmission limits
    # @constraint(
    #     ucmodel,
    #     TXCapTo[i = 1:ntrans, t = 1:ntimepoints],
    #     f[i, t] <= TFmax[i]
    # )
    # @constraint(
    #     ucmodel,
    #     TXCapFrom[i = 1:ntrans, t = 1:ntimepoints],
    #     f[i, t] >= -TFmax[i]
    # )

    # # DCOPF constraints
    # @constraint(
    #     ucmodel,
    #     DCOPTX[i = 1:ntrans, t = 1:ntimepoints],
    #     f[i, t] == TX[i] * (theta[i, t] - theta[i, t])
    # )

    # Power flow constraints
    # @constraint(
    #     ucmodel,
    #     PowerFlow[i = 1:ntrans, t = 1:ntimepoints],
    #     f[i, t] ==
    #     (transmap[i, 1] - transmap[i, 2]) *
    #     TX[i] *
    #     (guc[transmap[i, 1], t] - guc[transmap[i, 2], t])
    # )

    # Conventional generator capacity limits
    @constraint(
        ucmodel,
        UCCapU[i = 1:nucgen, t = 1:Horizon],
        guc[i, t] <= u[i, t] * GPmax[i]
    )
    @constraint(
        ucmodel,
        UCCapL[i = 1:nucgen, t = 1:Horizon],
        guc[i, t] >= u[i, t] * GPmin[i]
    )
    # Hydro and renewable capacity limits
    @constraint(
        ucmodel,
        HCap[i = 1:nhydro, t = 1:Horizon],
        gh[i, t] <= HAvail[i, t]
    )
    @constraint(
        ucmodel,
        ReCap[i = 1:nrenewable, t = 1:Horizon],
        gr[i, t] <= RAvail[i, t]
    )
    # Ramping limits
    @constraint(ucmodel, RUIni[i = 1:nucgen], guc[i, 1] - GPIni[i] <= GRU[i])
    @constraint(ucmodel, RDIni[i = 1:nucgen], GPIni[i] - guc[i, 1] <= GRD[i])
    @constraint(
        ucmodel,
        RU[i = 1:nucgen, t = 1:Horizon-1],
        guc[i, t+1] - guc[i, t] <= GRU[i]
    )
    @constraint(
        ucmodel,
        RD[i = 1:nucgen, t = 1:Horizon-1],
        guc[i, t] - guc[i, t+1] <= GRD[i]
    )

    # State transition constraints
    # @constraint(ucmodel, U0[i = 1:nucgen], u[i, 1] == u[i, ntimepoints])
    # Minimum up/down time constraints
    # @constraint(
    #     ucmodel,
    #     UT[i = 1:nucgen],
    #     u[i, 1] - U0[i] <= (1 - u[i, t]) * GUT[i]
    # )
    # if t > 1
    #     ### Minimum up/down time constraints
    #     for i in 1:nucgen
    #         @constraint(
    #             ucmodel,
    #             u[i, t] - u[i, t-1] <= (1 - u[i, t]) * GUT[i]
    #         )
    #         @constraint(
    #             ucmodel,
    #             u[i, t-1] - u[i, t] <= (1 - u[i, t-1]) * GDT[i]
    #         )
    #     end
    # end
    return ucmodel
end