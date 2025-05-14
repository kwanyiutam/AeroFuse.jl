module ConstraintDiagram

# Initialise packages used
using DataFrames
using Unitful
using ISAData


# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)


"""
    `point_perf` - A function to get the point performance of aircraft during climb, cruise, loiter or descend

    return the TW (or PW)
"""
function point_perf(;WS,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    # Mission information
    stage_row = findfirst(==("Stage"),df_mission[:,1])
    h_row = findfirst(==("Altitude"),df_mission[:,1])
    V∞_row = findfirst(==("Velocity"),df_mission[:,1])
    a_max_row = findfirst(==("Maximum Acceleration"),df_mission[:,1])
    gradient_row = findfirst(==("Gradient"),df_mission[:,1])
    engine_avail_row = findfirst(==("Available_Engines"),df_mission[:,1])
    α_row = findfirst(==("Cumulative Weight Fraction"),df_mission[:,1])
    β_row = findfirst(==("Engine Efficiency Scaling"),df_mission[:,1])
    n_row = findfirst(==("Load Factor"),df_mission[:,1])

    # Aircraft information
    η_prop = df_aircraft[findfirst(==("Engine Propeller Efficiency"),df_aircraft[:,1]),aircraft_idx]
    AR = df_aircraft[findfirst(==("Wing AR"),df_aircraft[:,1]),aircraft_idx]
    e = df_aircraft[findfirst(==("Oswald Efficiency"),df_aircraft[:,1]),aircraft_idx]
    CD0 = df_aircraft[findfirst(==("CD0"),df_aircraft[:,1]),aircraft_idx]
    engine_type = df_aircraft[findfirst(==("Engine Type"),df_aircraft[:,1]),aircraft_idx]
    number_of_engines = df_aircraft[findfirst(==("Number of Engines"),df_aircraft[:,1]),aircraft_idx]

    # Initialise a variable to save TW data
    TW = fill(0.0, length(WS), N_stages)

    idx = 1

    # For each column of mission
    for col in (ncol(df_mission)-N_stages+1):ncol(df_mission)
        stage = df_mission[stage_row,col]

        if stage in ["Climb","Cruise","Loiter","Descend"]
            # Get the needed information
            h = uconvert(u"m", df_mission[h_row,col])
            V∞ = uconvert(u"m/s", df_mission[V∞_row,col])
            a_max = uconvert(u"m/s^2", df_mission[a_max_row,col])
            engines_available = df_mission[engine_avail_row,col]
            gradient = df_mission[gradient_row,col]
            α = df_mission[α_row,col]
            β = df_mission[β_row,col]
            n = df_mission[n_row,col]
        
            ### Determine whether the input is power or thrust 
            if engine_type == "Jet"
                prop_constant = 1
            else
                prop_constant = V∞ / η_prop
            end

            ### Determine density
            (ρ,_,_,_) = ISAdata(h)

            # Special case for cruise and loiter, gradient would be zero
            if ismissing(gradient)
                gradient = 0.0
            end

            ### Terms in the point-performance
            g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
            C_climb = sin(atan(gradient / 100)) # Related to climb gradient
            C_speed = a_max / g # Related to change in speed
            C_zerol = (0.5 .* ρ .* (V∞ .^ 2) .* CD0) ./ (α .* WS)
            C_induc = (α .* n ^2 .* WS) ./ (0.5 .* ρ .* (V∞ .^ 2) .* pi .* AR .* e)

            TW[:,idx] = ustrip.((α ./ β) .* (number_of_engines ./ engines_available) .* prop_constant .* (C_climb .+ C_speed .+ C_zerol .+ C_induc))
        end
        idx += 1
    end

    return TW
end

"""
    `get_max_WS` - A function to get the maximum wing loading from the landing conditions

    return the TW
"""
function get_max_WS(;df_aircraft::DataFrame, aircraft_idx::Int,df_mission::DataFrame)
    row_idx = findfirst(==("Stage"),df_mission[:,1])
    landing_col = findall(==("Landing"),skipmissing(collect(df_mission[row_idx, :])))

    K_l = 1/df_aircraft[findfirst(==("Landing Distance Multiplier"), df_aircraft[:,1]),aircraft_idx] # Multiplier for landing distance (smaller)
    S_a = df_aircraft[findfirst(==("Glide Slope"), df_aircraft[:,1]), aircraft_idx] # Glide Slope
    CL_max = df_aircraft[findfirst(==("Wing CLmax"), df_aircraft[:,1]), aircraft_idx]
    ΔCL_max = df_aircraft[findfirst(==("Wing CLmax Landing Flaps"), df_aircraft[:,1]), aircraft_idx]
    CL_max_land = CL_max + ΔCL_max
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant

    WS_max = 999999
    (ρ_0,_,_,_) = ISAdata(0*u"m")

    for col in landing_col
        # Get landing distance
        LDA = df_mission[findfirst(==("Distance"), df_mission[:,1]),col] # Landing distance
        α = df_mission[findfirst(==("Max_Mass_Ratio"), df_mission[:,1]),col]
        reverser = df_mission[findfirst(==("Reverser"), df_mission[:,1]),col]
        (ρ,_,_,_) = ISAdata(df_mission[findfirst(==("Altitude"),df_mission[:,1]),col])
        σ = ρ / ρ_0

        if reverser == false
            K_r = 1
        else
            K_r = 0.66 # If reverser is allowed
        end

        WS = (((LDA * K_l - S_a)*(σ * CL_max_land)) ./ ((5/g) .* K_r .* α))

        WS_max = min(WS_max, ustrip(WS))
    end
    
    return WS_max*u"kg/m/s^2"
end

"""
    `quick_constraint` - A function to get the minimum T/W (or P/W for prop) required at a specific W/S (for quick mode)

    return the TW
"""
function quick_constraint(;WS_max,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    # Get the power/thrust from the conditions except for takeoff
    TW = point_perf(WS=[WS_max],df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages)

    # TODO: add one more for takeoff

    
    # Return the results from the search
    return maximum(TW)
end

end