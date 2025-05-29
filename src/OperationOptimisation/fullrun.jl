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


t1 = time()
(design_param,df_aircraft,N_aircraft,df_mission,N_stages,df_payload,N_payload,df_cost) = OperationOptimisation.design_init(aircraft_path="testing_ac.csv",payload_path="testing_payload.csv",mission_path="testing_mission.csv");

df_design = CSV.read("testing_sample.csv", DataFrame)
#df_design = OperationOptimisation.sampling_points(N_sample = 10, design_param = design_param, design_type = "Design")
#df_pert = OperationOptimisation.sampling_points(N_sample = 10, design_param = design_param, design_type = "Perturbations")
df_pert = CSV.read("testing_perturb.csv", DataFrame)
# By default it will return all cost elements, design geometry and optimised variables
data_save = ["CL_cruise","CDi_cruise","CDv_cruise","CDv_cruise_fuse","Block Time","LD_cruise","MTOW","Fuel Weight","Payload Weight","Empty Weight","Duration","Wing Span","Wing Area","HT Span","HT Area","HT Tail Arm","VT Span","VT Area","VT Tail Arm","Tmax","Fuel Tank Volume","Fuselage Length","Diameter"]

df_save = OperationOptimisation.design_start(design_param=design_param,df_design=df_design,df_pert=df_pert,df_aircraft=df_aircraft,N_aircraft=N_aircraft,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,df_cost=df_cost,data_save=data_save)


elapsed_time = time() - t1
print("Elasped Time: $elapsed_time seconds")
#print(df_save)

min_cost_idx = argmin(df_save[.!(isnan.(df_save[:,"Total Cost"])),"Total Cost"])

Plots.plot(
    dpi=500,
    aspect_ratio = :equal,
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

title!("Side View")

savefig("OptimalPlane_Side.png") 

Plots.plot(
    dpi=500,
    aspect_ratio = :equal,
    camera = (90, 0),
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

title!("Front View")

savefig("OptimalPlane_Front.png") 


Plots.plot(
    dpi=500,
    aspect_ratio =:equal,
    camera = (90, 90),
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

title!("Top View")


savefig("OptimalPlane_Top.png") 

detail_idx = findfirst(==("Cost Breakdown"),DataFrames.names(df_save))

save_data_design = Array(design_param[:,"Design Parameter"])

global keep_col = ["Total Cost","MTOW","Fuel Weight","Payload Weight","Empty Weight","CL_cruise","CDi_cruise","CDv_cruise","CDv_cruise_fuse","Block Time","LD_cruise","Wing Span","Wing Area","HT Span","HT Area","HT Tail Arm","VT Span","VT Area","VT Tail Arm","Tmax","Fuel Tank Volume","Fuselage Length","Diameter"]

keep_col = vcat(save_data_design,keep_col)

for row in 1:nrow(df_save)
    df_breakdown = df_save[row,detail_idx]

    for i in 1:nrow(df_breakdown)
        if startswith(df_breakdown[i, 1], "Aircraft Cost")
            df_breakdown[i, 1] = "Aircraft Cost"
        end
    end
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

try
    XLSX.writetable("TaperRatio.xlsx", df_save)
catch
    @warn "Failed to write table, creating something else!"
    XLSX.writetable("lbablbah1.xlsx", df_save)
end


