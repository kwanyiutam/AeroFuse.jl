using DataFrames
using Unitful
using CSV

# include the OperationOptimisation module
include("operationoptimisation.jl")

# include the CostModel module
include("costmodel.jl")

(design_param,df_aircraft,N_aircraft,df_mission,N_stages,df_payload,N_payload,df_cost) = OperationOptimisation.design_init();

df_design = CSV.read("./Reference/SampleDesign.csv", DataFrame)
#df_design = OperationOptimisation.sampling_points(N_sample = 10, design_param = design_param, design_type = "Design")
df_pert = OperationOptimisation.sampling_points(N_sample = 10, design_param = design_param, design_type = "Perturbations")

OperationOptimisation.design_start(design_param=design_param,df_design=df_design,df_pert=df_pert,df_aircraft=df_aircraft,N_aircraft=N_aircraft,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,df_cost=df_cost)

