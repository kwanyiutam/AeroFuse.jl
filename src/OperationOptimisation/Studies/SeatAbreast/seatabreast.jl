using DataFrames
using XLSX
using Plots
using Printf
using LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font, guidefont = font(14,plot_font), tickfont = font(12,plot_font), legendfont = font(8,plot_font), titlefont = font(16,plot_font), colorbar_titlefont = font(12,plot_font))

df_data = DataFrame(XLSX.readtable("./Studies/SeatAbreast/SeatAbreastInvestigate.xlsx",1))
#print(df_data)
cruise_range = 3500.0
pax = 250
after_body = 2.5
seat_abreast = 6
y = ["Cruise Range"]
x = ["Seat Abreast"]
metric = ["Total Cost"] # Total Cost Without Depreciation
flip_annotate = true
save = vcat(x,y,metric)

df_plot = df_data
#df_plot = df_data[df_data[:,"Cruise Range"]==cruise_range && df_data[:,"Number of Passengers"]==pax && df_data[:,"Afterbody Fineness Ratio"]==after_body,save]
#df_plot = df_data[df_data[:,"Seat Abreast"]==seat_abreast,save]
df_plot = filter(row -> (row["Number of Passengers"] == pax), df_data) # && row["Cruise Range"] == cruise_range
#&& df_data[:,"Seat Abreast"]==seat_abreast

df_plot_x = convert.(Int64,Array(df_plot[!,x]))
df_plot_y = convert.(Int64,Array(df_plot[!,y]))
df_plot_y_unique = unique(df_plot_y)
#df_plot_x = unique(convert.(Int64,Array(df_plot[!,x])))
#df_plot_y = unique(convert.(Int64,Array(df_plot[!,y])))
df_plot_metric = Array(df_plot[!,metric])
df_plot_metric[findall(==("NaN"),df_plot_metric)] .= NaN
df_plot_metric = convert.(Float64,df_plot_metric)
#df_plot_metric = reshape(df_plot_metric,(length(df_plot_y),length(df_plot_x)))
p = plot(dpi=500, legend = false)
metric_label = metric[1]
y_name = y[1]

for i in eachindex(df_plot_y_unique)
    y_i = df_plot_y_unique[i]
    idx = findall(==(y_i),df_plot_y)
    num = length(idx)
    x_i = df_plot_x[idx]
    z_i = df_plot_metric[idx]
    min_z_i = minimum(x for x in z_i if !isnan(x))
    z_i = (z_i .- min_z_i) ./ min_z_i * 100

    if flip_annotate == true
        offset = 0.1
        j = 1
        global success = true

        while j < num
            if isnan(z_i[j])
                j += 1
            else
                global x_annotate = x_i[j]
                global z_annotate = z_i[j]
                break
            end

            if j == num
                success = false
            end
        end
    else
        offset = -0.1
        j = num
        global success = true

        while j > 0
            if isnan(z_i[j])
                j -= 1
            else
                global x_annotate = x_i[j]
                global z_annotate = z_i[j]
                break
            end

            if j == 0
                success = false
            end
        end
    end

    if success == true
        plot!(p, x_i, z_i) #, linestyle = "Line&$y_i km&", show=true
        annotate!(p,x_annotate+offset,z_annotate,text("$y_i km",7))
    end
end
title!("$pax Passengers")
xlabel!("Seat Abreast")
ylabel!("Total Cost Using Full Aircraft Cost Model\n(% Difference From Minimum)") # Minus Depreciation

savefig("TotalCost_$pax.png")  # MinusDepreciation


### Plot for
p2 = plot(dpi=500,legend=:outerright)

x_plot = convert.(Int64, Array(df_data[:,"Number of Passengers"]))
y_plot = convert.(Float64, Array(df_data[:,"Cruise Range"]))
z_plot = convert.(Int64, Array(df_data[:,"Seat Abreast"]))
metric_plot = Array(df_data[:,metric])

# Get number of unique elements
z_plot_unique = unique(z_plot)
z_len = length(z_plot_unique)
z_max = maximum(z_plot_unique)
z_min = minimum(z_plot_unique)
total_len = length(z_plot)
itr = convert(Int64, total_len / z_len)

#my_colors = [cgrad(:phase, [0.01, 0.99])[z] for z ∈ range(0.0, 1.0, length = z_len)]
my_colors = get_color_palette(:default, plot_color(:white))

# Plot example markers for legend
for z_val in z_plot_unique
    scatter!(p2,
        [NaN], [NaN],  # dummy invisible point
        color = my_colors[z_val - z_min + 1],
        markersize = z_val,
        label = "$z_val-abreast"
    )
end
for i in 0:(itr-1)
    x_plot_i = x_plot[z_len*i+1]
    y_plot_i = y_plot[z_len*i+1]
    z_plot_i = z_plot[z_len*i+1:z_len*(i+1)]
    metric_plot_i = metric_plot[z_len*i+1:z_len*(i+1)]
    find_nan = findall(==("NaN"),metric_plot_i)
    metric_plot_i[find_nan] .= +Inf
    min_idx = argmin(skipmissing(metric_plot_i))
    z_plot_min = z_plot_i[min_idx]
    
    border = "black"

    if !isempty(find_nan)
        border = "red"
    end

    scatter!(p2, [x_plot_i], [y_plot_i], color = my_colors[z_plot_min-z_min+1],
        markersize = z_plot_min, label = "") #, markerstrokecolor = border, markerstrokewidth = 2
end

xlabel!("Number of Passengers")
ylabel!("Cruise Range (km)")
title!("Minimum Cost Configurations (Full Aircraft Model)") #"Minimum Cost Configurations (Without Depreciation)

#savefig("MinCostConfig.png")  # MinCostConfigNoDepreciate #MinCostConfigWithDepreciate #MinFuelConfig

p2