using DataFrames
using Unitful
using CSV
using XLSX
using Profile
using AeroFuse
using Plots
using LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font)

# include the OperationOptimisation module
include("operationoptimisation.jl")

# include the CostModel module
include("costmodel.jl")

(design_param,df_aircraft,N_aircraft,df_mission,N_stages,df_payload,N_payload,df_cost) = OperationOptimisation.design_init(aircraft_path="testing_ac.csv",payload_path="testing_payload.csv",mission_path="testing_mission.csv");

df_design = CSV.read("testing_sample.csv", DataFrame)
#df_design = OperationOptimisation.sampling_points(N_sample = 10, design_param = design_param, design_type = "Design")
df_pert = OperationOptimisation.sampling_points(N_sample = 10, design_param = design_param, design_type = "Perturbations")
#df_pert = CSV.read("testing_perturb.csv", DataFrame)
data_save = [""]

df_save = OperationOptimisation.design_start(design_param=design_param,df_design=df_design,df_pert=df_pert,df_aircraft=df_aircraft,N_aircraft=N_aircraft,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,df_cost=df_cost,data_save=data_save, mode = "Constraints")

print(df_save)

df_save = ustrip.(df_save)
#for col in 1:ncol(df_save)
#    df_save[:,col] = ustrip.(df_save[:,col])
#end

try
    XLSX.writetable("ConstraintDiagram.xlsx", df_save)
catch
    @warn "Failed to write table, creating something else!"
    XLSX.writetable("rsdifuhseuiprfghoweu.xlsx", df_save)
end

