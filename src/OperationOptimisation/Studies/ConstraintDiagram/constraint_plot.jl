using DataFrames
using XLSX
using Plots
using Printf
using LaTeXStrings
using Unitful

plot_font = "Computer Modern"
default(fontfamily=plot_font, guidefont = font(14,plot_font), tickfont = font(12,plot_font), legendfont = font(8,plot_font), titlefont = font(16,plot_font), colorbar_titlefont = font(12,plot_font))

df_data_min_fuel = DataFrame(XLSX.readtable("./Studies/ConstraintDiagram/ConstraintDiagram_MinFuel.xlsx",1))
df_data_min_cost = DataFrame(XLSX.readtable("./Studies/ConstraintDiagram/ConstraintDiagram_MinCost.xlsx",1))

x = ["WS"]
y = ["Cruise Velocity"]
metric = ["Total Cost"]
save_data = vcat(x,y,metric)

landing_distance = 1.5
filter!(row -> (row["Landing Distance"] == landing_distance), df_data_min_fuel) 
filter!(row -> (row["Landing Distance"] == landing_distance), df_data_min_cost) 
filter!(row -> (row["Cruise Velocity"] > 0.0), df_data_min_fuel) 
filter!(row -> (row["Cruise Velocity"] > 0.0), df_data_min_cost) 
filter!(row -> (row["WS"] > 440.0), df_data_min_fuel) 
filter!(row -> (row["WS"] > 440.0), df_data_min_cost) 

df_use = df_data_min_cost

p = plot(dpi=300,legend=:outerright)
plot!(p, df_use[:,"WS"], df_use[:,"TW Initial Climb"], label = "Initial Climb")
plot!(p, df_use[:,"WS"], df_data_min_cost[:,"TW Cruise to Destination"], linewidth = 3, label = "Cruise to Destination\n(Minimum Cost)")
plot!(p, df_use[:,"WS"], df_data_min_fuel[:,"TW Cruise to Destination"], linewidth = 3, label = "Cruise to Destination\n(Minimum Fuel)")
plot!(p, df_use[:,"WS"], df_use[:,"TW Descend to Destination"], label = "Descend to Destination")
plot!(p, df_use[:,"WS"], df_use[:,"TW Climb to Alternate"], label = "Climb to Alternate")
plot!(p, df_use[:,"WS"], df_use[:,"TW Cruise to Alternate"], label = "Cruise to Alternate")
plot!(p, df_use[:,"WS"], df_use[:,"TW Descend to Loiter"], label = "Descend to Loiter")
plot!(p, df_use[:,"WS"], df_use[:,"TW Loiter"], label = "Loiter")
plot!(p, df_use[:,"WS"], df_use[:,"TW Descend to Alternate"], label = "Descend to Alternate")
plot!(p, df_use[:,"WS"], df_use[:,"TW Ground Roll Takeoff"], label = "Ground Roll Takeoff")
plot!(p, df_use[:,"WS"], df_use[:,"TW Climb OEI Takeoff"], label = "Climb OEI Takeoff")
plot!(p, df_use[:,"WS"], df_use[:,"TW BFL Takeoff"], label = "BFL Takeoff")
vline!(p, df_use[:,"Landing WS_max"], label = "Landing")
xlabel!("Wing Loading W/S "*L"\mathrm{kg/m/s^2}")
ylabel!("Thrust to Weight Ratio T/W")
title!("Constraint Diagram")
#savefig("ConstraintDiagram.png")

p_cost = plot(dpi=300, ylim = (0, maximum(df_data_min_fuel[:,"Total Cost"])))
sNames_cost = ["Fuel Cost", "Flight Crew Cost", "Cabin Crew Cost", "Depreciation Cost"]
areaplot!(p_cost, df_data_min_cost[:,"WS"], [df_data_min_cost[:,sNames_cost[1]]' ; df_data_min_cost[:,sNames_cost[2]]' ; df_data_min_cost[:,sNames_cost[3]]' ; df_data_min_cost[:,sNames_cost[4]]']', labels=reshape(sNames_cost, (1,4)))
xlabel!(p_cost, "Wing Loading "*L"\mathrm{kg/m/s^2}")
ylabel!(p_cost, "Costs (USD)")
title!(p_cost, "Cost Distribution (Minimum Cost)")

#savefig("CostDistributionMinCost.png")

p_cost2 = plot(dpi=300, ylim = (0, maximum(df_data_min_fuel[:,"Total Cost"])))
sNames_cost = ["Fuel Cost", "Flight Crew Cost", "Cabin Crew Cost", "Depreciation Cost"]
areaplot!(p_cost2, df_data_min_fuel[:,"WS"], [df_data_min_fuel[:,sNames_cost[1]]' ; df_data_min_fuel[:,sNames_cost[2]]' ; df_data_min_fuel[:,sNames_cost[3]]' ; df_data_min_fuel[:,sNames_cost[4]]']', labels=reshape(sNames_cost, (1,4)))
xlabel!(p_cost2, "Wing Loading "*L"\mathrm{kg/m/s^2}")
ylabel!(p_cost2, "Costs (USD)")
title!(p_cost2, "Cost Distribution (Minimum Fuel)")

#savefig("CostDistributionMinFuel.png")

total_cost = df_data_min_cost[:,"Total Cost"]
p_cost3 = plot(dpi=300)
sNames_cost = ["Fuel Cost", "Flight Crew Cost", "Cabin Crew Cost", "Depreciation Cost"]
areaplot!(p_cost3, df_data_min_cost[:,"WS"], [(df_data_min_cost[:,sNames_cost[1]]./total_cost)' ; (df_data_min_cost[:,sNames_cost[2]]./total_cost)' ; (df_data_min_cost[:,sNames_cost[3]]./total_cost)' ; (df_data_min_cost[:,sNames_cost[4]]./total_cost)']', labels=reshape(sNames_cost, (1,4)))
xlabel!(p_cost3, "Wing Loading "*L"\mathrm{kg/m/s^2}")
ylabel!(p_cost3, "Costs Proportions")
title!(p_cost3, "Cost Distribution (Minimum Cost)")

#savefig("CostDistributionMinCost_prop.png")

total_cost = df_data_min_fuel[:,"Total Cost"]
p_cost4 = plot(dpi=300)
sNames_cost = ["Fuel Cost", "Flight Crew Cost", "Cabin Crew Cost", "Depreciation Cost"]
areaplot!(p_cost4, df_data_min_fuel[:,"WS"], [(df_data_min_fuel[:,sNames_cost[1]]./total_cost)' ; (df_data_min_fuel[:,sNames_cost[2]]./total_cost)' ; (df_data_min_fuel[:,sNames_cost[3]]./total_cost)' ; (df_data_min_fuel[:,sNames_cost[4]]./total_cost)']', labels=reshape(sNames_cost, (1,4)))
xlabel!(p_cost4, "Wing Loading "*L"\mathrm{kg/m/s^2}")
ylabel!(p_cost4, "Costs Proportions")
title!(p_cost4, "Cost Distribution (Minimum Fuel)")

#savefig("CostDistributionMinFuel_prop.png")


df_data_min_fuel = df_data_min_fuel[:,save_data]
df_data_min_cost = df_data_min_cost[:,save_data]

p = plot(dpi=300,top_margin = 5Plots.mm,right_marign = 2Plots.mm)

plot!(p, df_data_min_fuel[:,"WS"], df_data_min_fuel[:,"Cruise Velocity"], label="Minimum Fuel", linewidth = 2)
plot!(p, df_data_min_cost[:,"WS"], df_data_min_cost[:,"Cruise Velocity"], label="Minimum Cost", linewidth = 2)
xlabel!(p, "Wing Loading W/S "*L"\mathrm{kg/m/s^2}")
ylabel!(p, "Cruise Velocity knots")
title!(p, "Cruise Velocity Total Cost Comparisons")

cost_diff = (df_data_min_fuel[:,"Total Cost"] .- df_data_min_cost[:,"Total Cost"]) ./ df_data_min_fuel[:,"Total Cost"] .* 100
p2 = twinx()
bar_colors = [v ≥ 0 ? :green : :red for v in cost_diff]

bar!(p2, df_data_min_fuel[:,"WS"], cost_diff, color=bar_colors, linecolor=:transparent, legend = false)
ylabel!(p2, "Total Cost Percentage Improvement (%)")

savefig("VelocityCompare.png")  # MinCostConfigNoDepreciate #MinCostConfigWithDepreciate #MinFuelConfig
