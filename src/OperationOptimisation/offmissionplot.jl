using DataFrames
using XLSX
using CSV
using Plots
using Printf
using LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font, guidefont = font(14,plot_font), tickfont = font(12,plot_font), legendfont = font(8,plot_font), titlefont = font(16,plot_font), colorbar_titlefont = font(12,plot_font))

df_data = DataFrame(XLSX.readtable("./Studies/MultiMission/MissionDesignExplore.xlsx",1))
df_ops = CSV.read("testing_operational_design.csv", DataFrame)

#print(df_data)
cruise_range = 5000.0
pax = 150
after_body = 2.5
seat_abreast = 6
y = ["Number of Passengers"]
x = ["Cruise Range"]
metric = ["Average Flight Cost"]
save = vcat(x,y,metric)

df_plot = df_data
#df_plot = df_data[df_data[:,"Cruise Range"]==cruise_range && df_data[:,"Number of Passengers"]==pax && df_data[:,"Afterbody Fineness Ratio"]==after_body,save]

#&& df_data[:,"Seat Abreast"]==seat_abreast

df_plot_x = unique(convert.(Int64,Array(df_plot[!,x])))
df_plot_y = unique(convert.(Int64,Array(df_plot[!,y])))
df_plot_metric = Array(df_plot[!,metric])
df_plot_metric[findall(==("NaN"),df_plot_metric)] .= NaN
df_plot_metric = convert.(Float64,df_plot_metric)
df_plot_metric = reshape(df_plot_metric,(length(df_plot_y),length(df_plot_x)))

p = plot(dpi=500,legend=:outerbottom)

metric_label = metric[1]

contourf!(
    p,
    df_plot_x,
    df_plot_y,
    df_plot_metric,
    clabels=false,
    c = cgrad(:RdYlGn, rev = true),
    levels=20,
    alpha=0.7,
    left_margin = 20Plots.mm,
    right_margin = 15Plots.mm,
    colorbar_title= "\n\nAverage Flight Cost (USD)"
)

scatter!(
    p,
    df_ops[:,1],
    df_ops[:,2],
    label = "Off-Design Missions"
)

#x_opt = 14
#y_opt = 8
df_plot_metric[findall(isnan,df_plot_metric)] .= +Inf

opt_idx = argmin(df_plot_metric)
opt_range = df_plot_x[opt_idx[1]]
opt_pax = df_plot_y[opt_idx[2]]
println(opt_range)
println(opt_pax)
scatter!(p,[opt_range],[opt_pax], label = "Optimal Design Assumption: $opt_range km, $opt_pax passengers", markersize = 8, marker = :diamond)
#x_opt_idx = findfirst(==(x_opt), df_plot_x)
#y_opt_idx = findfirst(==(y_opt), df_plot_y)
#metric_opt = round(df_plot_metric[y_opt_idx,x_opt_idx],sigdigits = 3)
#scatter!(p,[x_opt],[y_opt],label="Selected Point ($x_opt,$y_opt), Error = $metric_opt%", markersize = 8)
xlabel!(p,x[1])
ylabel!(p,y[1])
title!(p,"Average Flight Cost for Different Design Assumptions")

p

savefig("MissionDesignExplore.png") 