using STESTS, JuMP, Gurobi, CSV, DataFrames

# Read data from .jld2 file
UCL, # hourly load for unit commitment, MW
genmap, # map of generators
GPmax, # maximum power output of generators
GPmin, # minimum power output of generators
GMC, # marginal cost of generators
GSMC, # segment marginal cost of generators
GINCPmax, # maximum power output of generator segments
transmap, # map of transmission lines
TX, # reactance of transmission lines
TFmax, # maximum power flow of transmission lines
GNLC, # no-load cost of generators
GRU, # ramp-up limit of generators
GRD, # ramp-down limit of generators
GSUC, # start-up cost of generators
GUT, # minimum up time of generators
GDT, # minimum down time of generators
GPini, # initial power output of generators
hydromap, # map of hydro units
HAvail, # availability of hydro units
renewablemap, # map of renewable units
RAvail,
storagemap,
EPC,
EPD,
Eeta,
ESOC,
ESOCini,
EDL,
EDHAvail,
EDRAvail = STESTS.read_jld2("./data/ADS2032.jld2")

UCHorizon = Int(24) # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
# NDay = Int(size(UCL, 1) / UCHorizon)
NDay = 1
EDSteps = Int(12) # number of 5-min intervals in a hour
# EDL = repeat(UCL, inner = (EDSteps, 1)) # load for economic dispatch, MW, repeat by ED steps w/o noise
#   # load for economic dispatch, MW, repeat by ED steps w/o noise
# EDRAvail = repeat(RAvail, inner = (EDSteps, 1)) # load for economic dispatch, MW, repeat by ED steps w/o noise
EDHorizon = Int(1) # optimization horizon for economic dispatch model, 1 without look-ahead, 12 with 1-hour look-ahead, 24 with 2-hour look-ahead, 48 with 4-hour look-ahead
EDGSMC = repeat(GSMC, outer = (1,1,EDHorizon)) # segment marginal cost of generators, repeat by UCHorizon
GSMC = repeat(GSMC, outer = (1,1,UCHorizon)) # segment marginal cost of generators, repeat by UCHorizon
#select first UCHorizon rows of UCL as initial input to unit commitment model
UCLInput = convert(Matrix{Float64}, UCL[1:UCHorizon, :]')
UCHAvailInput = convert(Matrix{Float64}, HAvail[1:UCHorizon, :]')
UCRAvailInput = convert(Matrix{Float64}, RAvail[1:UCHorizon, :]')
#select first EDHorizon rows of EDL as initial input to economic dispatch model
EDLInput = convert(Matrix{Float64}, EDL[1:EDHorizon, :]')
EDHAvailInput = convert(Matrix{Float64}, EDHAvail[1:EDHorizon, :]')
EDRAvailInput = convert(Matrix{Float64}, EDRAvail[1:EDHorizon, :]')
UInput = convert(Array{Int64,1}, GPini .!= 0) # initial status of generators, 1 if on, 0 if off

# Formulate unit commitment model
ucmodel = STESTS.unitcommitment(
    UCLInput,
    genmap,
    GPmax,
    GPmin,
    GMC,
    GSMC,
    GINCPmax,
    transmap,
    TX,
    TFmax,
    GNLC,
    GRU,
    GRD,
    GSUC,
    GUT,
    GDT,
    GPini,
    UInput,
    hydromap,
    UCHAvailInput,
    renewablemap,
    UCRAvailInput,
    storagemap,
    EPC,
    EPD,
    Eeta,
    ESOC,
    ESOCini,
    Horizon = UCHorizon, # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
    VOLL = 1000.0, # value of lost load, $/MWh
    RM = 0.03, # reserve margin, 6% of peak load
)

# Edit unit commitment model here
# set optimizer, set add_bridges = false if model is supported by solver
set_optimizer(ucmodel, Gurobi.Optimizer, add_bridges = false)
# # modify objective function
# @objective(ucmodel, Min, 0.0)
# # modify or add constraints
# @constraint(ucmodel, 0.0 <= ucmodel[:P][1,1] <= 0.0)

#  Formulate economic dispatch model
edmodel = STESTS.economicdispatch(
    EDLInput,
    genmap,
    GPmax,
    GPmin,
    GMC,
    EDGSMC,
    GINCPmax,
    transmap,
    TX,
    TFmax,
    GNLC,
    GRU,
    GRD,
    GPini,
    hydromap,
    EDHAvailInput,
    renewablemap,
    EDRAvailInput,
    UInput,
    storagemap,
    EPC,
    EPD,
    Eeta,
    ESOC,
    ESOCini,
    Horizon = EDHorizon,
    Steps = EDSteps, # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
    VOLL = 1000.0, # value of lost load, $/MWh
)

# Edit economic dispatch model here
# set optimizer, set add_bridges = false if model is supported by solver
set_optimizer(edmodel, Gurobi.Optimizer, add_bridges = false)
# # modify objective function
# @objective(edmodel, Min, 0.0)
# # modify or add constraints
# @constraint(edmodel, 0.0 <= edmodel[:P][1,1] <= 0.0)

# Solve
timesolve = @elapsed begin
    UCcost, UCnetgen, UCgen, EDcost = STESTS.solving(
        NDay,
        UCL,
        HAvail,
        RAvail,
        GRU,
        GRD,
        GPmax,
        GPmin,
        GPini,
        GMC,
        GSMC,
        GNLC,
        GSUC,
        EDL,
        EDHAvail,
        EDRAvail,
        ucmodel,
        edmodel,
        UCHorizon = UCHorizon,
        EDHorizon = EDHorizon,
        EDSteps = EDSteps,
        VOLL = 1000.0,
    )
end
@info "Solving took $timesolve seconds."

CSV.write("output/UCgentotal.csv", DataFrame(UCgen, :auto))
CSV.write("output/UCnetgentotal.csv", DataFrame(UCnetgen, :auto))

println("The UC cost is: ", sum(UCcost))
println("The ED cost is: ", sum(EDcost))
