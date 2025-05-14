using DataFrames
using Unitful
using CSV

# include the OperationOptimisation module
include("operationoptimisation.jl")

# include the CostModel module
include("costmodel.jl")

(design_param,df_aircraft,N_aircraft,df_mission,N_stages,df_payload,N_payload,df_cost) = OperationOptimisation.design_init(aircraft_path = "bogus.csv");
df_design = CSV.read("./Reference/SampleDesign.csv", DataFrame)

@profview OperationOptimisation.design_start(design_param=design_param,df_design=df_design,df_aircraft=df_aircraft,N_aircraft=N_aircraft,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,df_cost=df_cost)

