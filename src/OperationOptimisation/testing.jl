using DataFrames
using Unitful
using CSV
using XLSX
using Profile
using AeroFuse
using Plots

# include the OperationOptimisation module
include("operationoptimisation.jl")

# include the CostModel module
include("costmodel.jl")

(design_param,df_aircraft,N_aircraft,df_mission,N_stages,df_payload,N_payload,df_cost) = OperationOptimisation.design_init(aircraft_path="testing_ac.csv",payload_path="testing_payload.csv",mission_path="testing_mission.csv");

df_design = CSV.read("testing_sample.csv", DataFrame)
#df_design = OperationOptimisation.sampling_points(N_sample = 10, design_param = design_param, design_type = "Design")
#df_pert = OperationOptimisation.sampling_points(N_sample = 10, design_param = design_param, design_type = "Perturbations")
df_pert = CSV.read("testing_perturb.csv", DataFrame)
# By default it will return all cost elements, design geometry and optimised variables
data_save = ["CL_cruise","CDi_cruise","CDv_cruise","LD_cruise","MTOW","Fuel Weight","Payload Weight","Empty Weight","Duration","Wing Span","Wing Area","HT Span","HT Area","HT Tail Arm","VT Span","VT Area","VT Tail Arm","Tmax","Fuel Tank Volume","Fuselage Length","Diameter"]

df_save = OperationOptimisation.design_start(design_param=design_param,df_design=df_design,df_pert=df_pert,df_aircraft=df_aircraft,N_aircraft=N_aircraft,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,df_cost=df_cost,data_save=data_save)

#print(df_save)

min_cost_idx = argmin(df_save[.!(isnan.(df_save[:,"Total Cost"])),"Total Cost"])

Plots.plot(
    aspect_ratio = 1,
    camera = (0, 0),
    zlim = span(df_save[min_cost_idx,"Wing Mesh"]) .* (-0.5, 0.5),
    size = (800, 600)
)
Plots.plot!(df_save[min_cost_idx,"Wing Mesh"], label = "Wing")
Plots.plot!(df_save[min_cost_idx,"HT Mesh"], label = "HT")
Plots.plot!(df_save[min_cost_idx,"VT Mesh"], label = "VT")
Plots.plot!(df_save[min_cost_idx,"Fuselage Shape"], label = "Fuselage")

for i in eachindex(df_save[min_cost_idx,"Engine Shape"])
    Plots.plot!(df_save[min_cost_idx,"Engine Shape"][i], label = "Engine $i")
end

savefig("OptimalPlane.png") 

detail_idx = findfirst(==("Cost Breakdown"),DataFrames.names(df_save))

global keep_col = ["Spanwise Discretisation","Chordwise Discretisation","Total Cost","MTOW","Fuel Weight","Payload Weight","Empty Weight","CL_cruise","CDi_cruise","CDv_cruise","LD_cruise","Wing Span","Wing Area","HT Span","HT Area","HT Tail Arm","VT Span","VT Area","VT Tail Arm","Tmax","Fuel Tank Volume","Fuselage Length","Diameter"]

for row in 1:nrow(df_save)
    df_breakdown = df_save[row,detail_idx]
    cost_component = df_breakdown[:,1]

    # Initialise columns
    if row == 1
        global keep_col = vcat(keep_col, cost_component)
        for i in 1:nrow(df_breakdown)
            colname = cost_component[i]
            df_save[!,colname] = fill(0.0, nrow(df_save))
        end
    end

    cost_detail = df_breakdown[:,"Cost (USD)"]

    for i in 1:nrow(df_breakdown)
        df_save[row,cost_component[i]] = cost_detail[i]
    end
end

df_save = df_save[:,keep_col]

for row in 1:nrow(df_save)
    for col in 1:ncol(df_save)
        df_save[row,col] = ustrip.(df_save[row,col])
    end
end

XLSX.writetable("TestingOutput_again.xlsx", df_save)


