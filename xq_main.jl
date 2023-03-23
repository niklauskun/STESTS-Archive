using Pkg
Pkg.instantiate()
using JuMP
using Gurobi
using MAT
using Printf
using DataFrames
using CSV
# using Plots
using JLD
using DelimitedFiles
using Statistics, StatsBase
using Random, Distributions
include("xq_fcn.jl")


# load the data fileln
filename = "DA1108_1seg_MP" # price zone name
fileln = matopen(string("./Data/", filename, ".mat"))
sys_pa = read(fileln, "sys_para")
close(fileln)


############################################################
## Simulation setting
############################################################
# System setting
detT = 1; # Time between two scheduling periods [hour]
T = 24; # Time periods. (Total time = T*detT)
Tw = 4; # Time window for single-slot dispatch. (Window time = Tw * detT)

Wc = [0.5] * 1.3e4; # wind capacity [times]
Wu = [20] / 100; # real-time wind uncertainty [%]
ES_Pc = [1 5 10 20 50 80 100 200]; # stoarge capacity [MW]
E_Pu = [20]; # the deviation (assume Guassian) used to consider bid uncertainty.

num_sc = 5; # number of total wind/demand scenarios
num_Wc = length(Wc); # number of wind capacity scenarios [MW]
num_Wu = length(Wu); # number of wind uncertainty scenarios [%]
num_Ec = length(ES_Pc); # number of storage capacity [MW]
num_Eu = length(E_Pu); # number of storage forecast of price uncertainty [$/MWh]

# Storage setting
Pr = 0.25; # power rating
Pmax = 1; # power 
ed = 0.01;   # SoC sample granularity
Ne = Int64(floor(1 / ed) + 1); # number of SOC samples
ef = 0.0;    # final SoC target level, use 0 if none
num_Eseg = 1;# value function segment number
e0 = 0'; # initial e0 [p.u.]
eta = 0.9; # storage efficiency [p.u.]
c_MC = 20; # marginal discharge cost [$/MWh]
Ep = 50;      # penalty cost if ES cannot finish, $/MWh

############################################################
## transform data
############################################################
iS = 1;
iWc = 1;
iWu = 1;
iPc = 1;
iPu = 1;
# load parameters

E_Ec = ES_Pc ./ Pr;
Pd_max = Pmax * ES_Pc[iPc];
Pc_max = Pmax * ES_Pc[iPc];
numG = length(sys_pa["g_max"]);
W = sys_pa["wind_after_kmeans"][:,iS] .* Wc[iWc];
D = sys_pa["load_after_kmeans"][:,iS];
Tup = sys_pa["gen"][:,4];
Tdn = sys_pa["gen"][:,5];
G_min = sys_pa["gen"][:,2];
G_max = sys_pa["gen"][:,1];
R_max = sys_pa["gen"][:,6];
c1 = sys_pa["gen"][:,9];
c2 = sys_pa["gen"][:,10];
c_NLC = sys_pa["gen"][:,8];
c_SC = sys_pa["gen"][:,7];
tao_up = zeros(Int64, numG, T);
tao_dn = zeros(Int64, numG, T);
for t = 1:T
  for s = 1:numG
    tao_up[s,t] = Int64(round(max(t - Tup[s] +1,1)));
    tao_dn[s,t] = Int64(round(max(t - Tup[s] +1,1)));
  end
end
############################################################
## Block 1: Unit commitment
############################################################
for test = 1:1
  ############################################################
  ## UC problem
  ############################################################
  # initialize optimization model
  model = Model(Gurobi.Optimizer)
  set_silent(model) # no outputs
  # Storage: discharge and charge power, SoC level
  @variable(model, p_d[1:T], lower_bound = 0) # discharge power [MW]
  @variable(model, p_c[1:T], lower_bound = 0) # charge power [MW]
  @variable(model, e[1:T], lower_bound = 0)   # SoC [MWh]
  # Generator: generation and unit commitment
  @variable(model, g[1:numG, 1:T], lower_bound = 0) # generator power [MW]
  @variable(model, u[1:numG, 1:T], Bin) # logic variable, operation status
  @variable(model, v[1:numG, 1:T], Bin) # logic variable, strat-up status
  @variable(model, z[1:numG, 1:T], Bin) # logic variable
  # Wind: 
  @variable(model, w[1:T], lower_bound = 0) # wind generation [MW]
  # Reserve:
  @variable(model, r[1:T], lower_bound = 0) # system reserve [MW]
  #######################################################
  ####         Constraints
  #######################################################
  # generation max and min limits
  @constraint(model, C_Gmax[t=1:T], g[:,t] .<= G_max .* u[:,t] )
  @constraint(model, C_Gmin[t=1:T], g[:,t] .>= G_min .* u[:,t] )
  # No ramping limits here
  # UC start-up / shut-down logic constraints
  @constraint(model, C_Logic1[s=1:numG, t=2:T], v[s,t] - z[s,t] == u[s,t] - u[s,t-1] )
  @constraint(model, C_Logic2[s=1:numG, t=2:T], v[s,t] - z[s,t] <= 1 )
  # Minimum up time constraint
  @constraint(model, C_Minup[s=1:numG,t=1:T], sum(v[s,tao_up[s,t]:t]) <= u[s,t]);
  @constraint(model, C_Mindn[s=1:numG,t=1:T], sum(v[s,tao_dn[s,t]:t]) <= 1-u[s,t]);
  # system power balance constraint
  @constraint(model, C_bal[t=1:T], sum(g[:,t]) + w[t] + sum(p_d[t]) == D[t] + sum(p_c[t]));
  # system wind generation constraint
  @constraint(model, C_w[t=1:T], w[t] <= W[t]);
  # storage maximum power for discharging/charging
  @constraint(model, C_PdMax[t=1:T], p_d[t]./Pd_max <= 1 )
  @constraint(model, C_PcMax[t=1:T], p_c[t]./Pc_max <= 1 )
  # SoC evolution
  @constraint(model, C_SoC1[t=2:T], e[t] - e[t-1] == detT * (p_c[t] *eta - p_d[t]/eta)) ;
  @constraint(model, C_SoC2,        e[1] - e0 == detT * (c_MC*p_c[1]* eta - p_d[1]/eta)) ;
  # System reserve
  @constraint(model, C_Res1[t=1:T], r[t] >= 0.2*w[t]) ;
  @constraint(model, C_Res2[t=1:T], r[t] <= sum(G_max .* u[:,t] - g[:,t])) ; 
  @constraint(model, C_Res3[t=1:T], r[t] <= sum(R_max .* u[:,t])) ;
  # Objective function
  @objective(model, Min, sum(sum(c1[s] * g[s,t] + c2[s] * (g[s,t]^2) 
        + u[s,t] * c_NLC[s] +  v[s,t] * c_SC[s] for s = 1:numG) for t = 1:T)) 
  optimize!(model)
  ############################################################
  ## Print results
  ############################################################
  R_Gcost = value(sum(sum(c1[s] * g[s,t] + c2[s] * (g[s,t]^2) for s = 1:numG) for t = 1:T)) # total generation costs
  R_Ecost = value(sum(p_d .* detT .* c_MC))  # total storage degradation costs
  R_Genergy = value(sum(sum(g))) # total generator energy
  R_EDenergy = value(sum(p_d))  # total storage discharge energy
  R_ECenergy = value(sum(p_c))  # total storage charge energy
  R_Wenergy = value(sum(w))  # total wind generation energy
  R_Lenergy = value(sum(D))  # total load demand
  # save integer variables
  value_u = value.(u);
  value_v = value.(v);
  value_z = value.(z);
  @printf("Generator cost %d, Storage degradation cost %d \n", R_Gcost, R_Ecost)
  @printf("generation energy %d [MWh], storage discharge energy %d [MWh], wind generation %d [MWh] \n", R_Genergy, R_EDenergy,R_Wenergy)
  @printf("Load energy %d [MWh], imbalance %d [MWh] \n", R_Lenergy, R_Genergy+R_EDenergy+R_Wenergy-R_Lenergy-R_ECenergy)
  @printf("OptStatus: %s \n", termination_status(model))


  ############################################################
  ## Dual problem
  ############################################################
  # initialize optimization model
  model2 = Model(Gurobi.Optimizer)
  set_silent(model2) # no outputs
  # discharge and charge power, SoC level
  @variable(model2, p_d[1:T], lower_bound = 0) # discharge power [MW]
  @variable(model2, p_c[1:T], lower_bound = 0) # charge power [MW]
  @variable(model2, e[1:T], lower_bound = 0)   # SoC [MWh]
  # energy 
  @variable(model2, g[1:numG, 1:T], lower_bound = 0) # generator power [MW]
  # wind variable
  @variable(model2, w[1:T], lower_bound = 0) # wind generation [MW]
  # system variable
  @variable(model2, r[1:T], lower_bound = 0) # system reserve [MW]
  #######################################################
  ####         Constraints
  #######################################################
  # generation max and min limits
  @constraint(model2, C_Gmax[t=1:T], g[:,t] .<= G_max .* value_u[:,t] )
  @constraint(model2, C_Gmin[t=1:T], g[:,t] .>= G_min .* value_u[:,t] )
  # We do not have ramping limits here

  # system power balance constraint
  @constraint(model2, C_bal[t=1:T], sum(g[:,t]) + w[t] + sum(p_d[t]) == D[t] + sum(p_c[t]));
  # system wind generation constraint
  @constraint(model2, C_w[t=1:T], w[t] <= W[t]);
  # storage maximum power for discharging/charging
  @constraint(model2, C_PdMax[t=1:T], p_d[t]./Pd_max <= 1 )
  @constraint(model2, C_PcMax[t=1:T], p_c[t]./Pc_max <= 1 )
  # SoC evolution
  @constraint(model2, C_SoC1[t=2:T], e[t] - e[t-1] == detT * (p_c[t] *eta - p_d[t]/eta)) ;
  @constraint(model2, C_SoC2,        e[1] - e0 == detT * (c_MC*p_c[1]* eta - p_d[1]/eta)) ;
  # System reserve
  @constraint(model2, C_Res1[t=1:T], r[t] >= 0.2*w[t]) ;
  @constraint(model2, C_Res2[t=1:T], r[t] <= sum(G_max .* value_u[:,t] - g[:,t])) ; 
  @constraint(model2, C_Res3[t=1:T], r[t] <= sum(R_max .* value_u[:,t])) ; 
  # Objective function
  @objective(model2, Min, sum(sum(c1[s] * g[s,t] + c2[s] * (g[s,t]^2) + 
          value_u[s,t] * c_NLC[s] +  value_v[s,t] * c_SC[s] for s = 1:numG) for t = 1:T)) 
  optimize!(model2)
  LMP = zeros(T,1);
  if dual_status(model2) == FEASIBLE_POINT
    LMP = dual.(C_bal);
    println("  dual solution: C_bal = ", dual.(C_bal))
  else
    println("  No dual solution")
end



############################################################
## Block 2: Value function calculation
############################################################
v_discharge, v_charge =  value_fcn_calcu(LMP ,num_Eseg, Ne, T, c_MC, ES_Pc[iPc], eta, ed, ef, E_Pu[iPu])

println("  discharge bid = ", v_discharge)
println("  charge bid = ", v_charge)
end