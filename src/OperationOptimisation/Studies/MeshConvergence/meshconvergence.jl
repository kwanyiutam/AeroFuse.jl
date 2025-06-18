using DataFrames
using XLSX
using Plots
using Printf
using LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font, guidefont = font(14,plot_font), tickfont = font(12,plot_font), legendfont = font(8,plot_font), titlefont = font(16,plot_font), colorbar_titlefont = font(12,plot_font))

df_data = DataFrame(XLSX.readtable("./Studies/MeshConvergence/ConvergenceStudy.xlsx",1))
#print(df_data)
cruise_range = 5000.0
pax = 150
after_body = 2.5
seat_abreast = 6
y = ["Spanwise Discretisation"]
x = ["Chordwise Discretisation"]
metric = ["CDi_cruise"]
save = vcat(x,y,metric)

df_plot = df_data
#df_plot = df_data[df_data[:,"Cruise Range"]==cruise_range && df_data[:,"Number of Passengers"]==pax && df_data[:,"Afterbody Fineness Ratio"]==after_body,save]

#&& df_data[:,"Seat Abreast"]==seat_abreast

df_plot_x = unique(convert.(Int64,Array(df_plot[!,x])))
df_plot_y = unique(convert.(Int64,Array(df_plot[!,y])))
df_plot_metric = Array(df_plot[!,metric])
df_plot_metric[findall(==("NaN"),df_plot_metric)] .= NaN
df_plot_metric = convert.(Float64,df_plot_metric)
ref_metric = df_plot_metric[length(df_plot_metric)]
df_plot_metric = (df_plot_metric .- ref_metric) ./ ref_metric * 100 # Percentage difference
df_plot_metric = reshape(df_plot_metric,(length(df_plot_y),length(df_plot_x)))
p = plot(dpi=500, legend = :outerbottom)

metric_label = metric[1]

contourf!(
    df_plot_x,
    df_plot_y,
    df_plot_metric,
    clabels=false,
    color = get_color_palette(:speed, plot_color(:white)),
    levels=20,
    alpha=0.7,
    right_margin = 5Plots.mm,
    colorbar_title= "\n % Difference from Most Refined "*L"\mathrm{C_{Di,cruise}}"
)

x_opt = 14
y_opt = 8
x_opt_idx = findfirst(==(x_opt), df_plot_x)
y_opt_idx = findfirst(==(y_opt), df_plot_y)
metric_opt = round(df_plot_metric[y_opt_idx,x_opt_idx],sigdigits = 3)
scatter!(p,[x_opt],[y_opt],label="Selected Point ($x_opt,$y_opt), Error = $metric_opt%", markersize = 8)
xlabel!(p,x[1])
ylabel!(p,y[1])
title!(p,"Mesh Convergence Study")

p

savefig("MeshConvergence.png") 