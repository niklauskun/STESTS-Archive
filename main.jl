using STESTS, JuMP, Gurobi, CSV, DataFrames, Statistics

# Read data from .jld2 file 
params = STESTS.read_jld2("./data/ADS2032_7RegionNoise_4hrBES_5GWBES.jld2")
model_filenames =
    ["models/WEST_1.jld2", "models/WEST_2.jld2", "models/WEST_3.jld2"]

strategic = false
RM = 0.03
VOLL = 9000.0
# PriceCap = 1000.0
UCHorizon = Int(25) # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
EDHorizon = Int(13) # optimization horizon for economic dispatch model, 1 without look-ahead, 12 with 1-hour look-ahead
NDay = 364
EDSteps = Int(12) # number of 5-min intervals in a hour
ESSeg = Int(1)
PriceCap = repeat(
    repeat((range(220, stop = 1000, length = 40))', outer = (7, 1)),
    outer = (1, 1, EDHorizon),
)
FuelAdjustment = 1.5
ErrorAdjustment = 0.25
LoadAdjustment = 1.0

output_folder =
    "output/DecUpdate/UC" *
    "$UCHorizon" *
    "ED" *
    "$EDHorizon" *
    "_Strategic_" *
    "$strategic" *
    "_Seg" *
    "$ESSeg" *
    "_Load" *
    "$LoadAdjustment" *
    "_Fuel" *
    "$FuelAdjustment" *
    "_Error" *
    "$ErrorAdjustment" *
    "_5GWBES_1yr_emergency"
mkpath(output_folder)

DADBidsSingle = [
    150,
    150,
    145,
    140,
    135,
    130,
    125,
    140,
    150,
    150,
    150,
    150,
    150,
    145,
    140,
    120,
    100,
    80,
    60,
    50,
    80,
    110,
    140,
    150,
    150,
    150,
    145,
    140,
    135,
    130,
    125,
    140,
    150,
    150,
    150,
    150,
    150,
    145,
    140,
    120,
    100,
    80,
    60,
    50,
    80,
    110,
    140,
    150,
    150,
    150,
    145,
    140,
    135,
    130,
    125,
    140,
    150,
    150,
    150,
    150,
    150,
    145,
    140,
    120,
    100,
    80,
    60,
    50,
    80,
    110,
    140,
    150,
]
DACBidsSingle = [
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -35,
    -20,
    -5,
    10,
    25,
    10,
    -5,
    -20,
    -35,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -35,
    -20,
    -5,
    10,
    25,
    10,
    -5,
    -20,
    -35,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -35,
    -20,
    -5,
    10,
    25,
    10,
    -5,
    -20,
    -35,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
    -50,
]
DADBids = repeat(DADBidsSingle', size(params.storagemap, 1), 1)
DACBids = repeat(DACBidsSingle', size(params.storagemap, 1), 1)
RTDBids = repeat(DADBids, inner = (1, EDSteps))
RTCBids = repeat(DACBids, inner = (1, EDSteps))

bidmodels = STESTS.loadbidmodels(model_filenames)
storagebidmodels =
    STESTS.assign_models_to_storages(bidmodels, size(params.storagemap, 1))

# Formulate unit commitment model
ucmodel = STESTS.unitcommitment(
    params,
    Horizon = UCHorizon, # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
    VOLL = VOLL, # value of lost load, $/MWh
    RM = RM, # reserve margin
    FuelAdjustment = FuelAdjustment,
)

# Edit unit commitment model here
# set optimizer, set add_bridges = false if model is supported by solver
set_optimizer(ucmodel, Gurobi.Optimizer, add_bridges = false)
set_optimizer_attribute(ucmodel, "OutputFlag", 0)
# # modify objective function
# @objective(ucmodel, Min, 0.0)
# # modify or add constraints
# @constraint(ucmodel, 0.0 <= ucmodel[:P][1,1] <= 0.0)

ucpmodel = STESTS.unitcommitmentprice(
    params,
    Horizon = UCHorizon, # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
    VOLL = VOLL, # value of lost load, $/MWh
    RM = RM, # reserve margin
    FuelAdjustment = FuelAdjustment,
)

# Edit unit commitment model here
# set optimizer, set add_bridges = false if model is supported by solver
set_optimizer(ucpmodel, Gurobi.Optimizer, add_bridges = false)
set_optimizer_attribute(ucpmodel, "OutputFlag", 0)
# # modify objective function
# @objective(ucpmodel, Min, 0.0)
# # modify or add constraints
# @constraint(ucpmodel, 0.0 <= ucpmodel[:P][1,1] <= 0.0)

#  Formulate economic dispatch model
edmodel = STESTS.economicdispatch(
    params,
    PriceCap, # value of lost load, $/MWh
    ESSeg = ESSeg,
    Horizon = EDHorizon,
    Steps = EDSteps, # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
    FuelAdjustment = FuelAdjustment,
)

# Edit economic dispatch model here
# set optimizer, set add_bridges = false if model is supported by solver
set_optimizer(edmodel, Gurobi.Optimizer, add_bridges = false)
set_optimizer_attribute(edmodel, "OutputFlag", 0)
# # modify objective function
# @objective(edmodel, Min, 0.0)
# # modify or add constraints
# @constraint(edmodel, 0.0 <= edmodel[:P][1,1] <= 0.0)

# Solve
timesolve = @elapsed begin
    UCcost, EDcost = STESTS.solving(
        params,
        NDay,
        strategic,
        DADBids,
        DACBids,
        RTDBids,
        RTCBids,
        ucmodel,
        ucpmodel,
        edmodel,
        storagebidmodels,
        output_folder,
        PriceCap,
        ESSeg = ESSeg,
        UCHorizon = UCHorizon,
        EDHorizon = EDHorizon,
        EDSteps = EDSteps,
        VOLL = VOLL,
        RM = RM,
        ErrorAdjustment = ErrorAdjustment,
        LoadAdjustment = LoadAdjustment,
    )
end
@info "Solving took $timesolve seconds."

println("The UC cost is: ", sum(UCcost))
println("The ED cost is: ", sum(EDcost))
