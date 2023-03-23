using STESTS, JuMP, Gurobi

# Read data from .jld2 file
D,
genmap,
GPmax,
GPmin,
GMC,
transmap,
TX,
TFmax,
GNLC,
GRU,
GRD,
GSUC,
GUT,
GDT,
GPIni,
hydromap,
HAvail,
renewablemap,
RAvail = STESTS.read_jld2("./data/WECC240.jld2")

UCHorizon = Int(24) # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
NDay = Int(size(D, 1) / UCHorizon)
#select first HCHorizon rows of D as initial input to unit commitment model
DInput = convert(Matrix{Float64}, D[1:UCHorizon, :]')

# Formulate unit commitment model
ucmodel = STESTS.unitcommitment(
    DInput,
    genmap,
    GPmax,
    GPmin,
    GMC,
    transmap,
    TX,
    TFmax,
    GNLC,
    GRU,
    GRD,
    GSUC,
    GUT,
    GDT,
    GPIni,
    hydromap,
    HAvail,
    renewablemap,
    RAvail,
    Horizon = UCHorizon, # optimization horizon for unit commitment model, 24 hours for WECC data, 4 hours for 3-bus test data
    VOLL = 9000.0, # value of lost load, $/MWh
)

# Edit unit commitment model here
# set optimizer, set add_bridges = false if model is supported by solver
set_optimizer(ucmodel, Gurobi.Optimizer, add_bridges = false)
# # modify objective function
# @objective(ucmodel, Min, 0.0)
# # modify or add constraints
# @constraint(ucmodel, 0.0 <= ucmodel[:P][1,1] <= 0.0)
# return latex formulation of ucmodel
latex_formulation(ucmodel)

# Formulate economic dispatch model
# edmodel = STESTS.economicdispatch(
#     c = data[1],
#     a = data[2],
#     p_max = data[3],
#     p_min = data[4],
#     demand = data[5],
# )

# Edit economic dispatch model here

# Solve
cost = STESTS.solving(1, D, ucmodel, Horizon = UCHorizon)

# optimize!(ucmodel)

# # Check status
# if termination_status(ucmodel) != MOI.OPTIMAL
#     error("Optimization failed.")
# end
# optimize!(edmodel)