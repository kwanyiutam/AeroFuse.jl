module ConstraintDiagram

# Initialise packages used
using DataFrames
using Unitful
using ISAData

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the ParseData module
include("parsedata.jl")

# Include the InputValidate module
include("inputvalidate.jl")

"""
    `point_perf` - A function to get the point performance of aircraft during climb, cruise, loiter or descend

    return the TW (or PW)
"""
function point_perf(;df_WS,velocity_list::DataFrame,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    # Mission information
    stage_row = findfirst(==("Stage"),df_mission[:,1])
    h_row = findfirst(==("Altitude"),df_mission[:,1])
    V∞_row = findfirst(==("Velocity"),df_mission[:,1])
    a_max_row = findfirst(==("Maximum Acceleration"),df_mission[:,1])
    gradient_row = findfirst(==("Gradient"),df_mission[:,1])
    engine_avail_row = findfirst(==("Available_Engines"),df_mission[:,1])
    α_row = findfirst(==("α"),df_mission[:,1])
    β_row = findfirst(==("Engine Efficiency Scaling"),df_mission[:,1])
    n_row = findfirst(==("Load Factor"),df_mission[:,1])
    WS = df_WS[:,1]

    # These are the columns where the velocity should be taken from WS instead!
    key_idx = velocity_list[:,"Saved Column"]

    # Aircraft information
    η_prop = InputValidate.get_value(df_aircraft,"Engine Propeller Efficiency",aircraft_idx)
    AR = InputValidate.get_value(df_aircraft,"Wing AR",aircraft_idx)
    e = InputValidate.get_value(df_aircraft,"Oswald Efficiency",aircraft_idx)
    CD0 = InputValidate.get_value(df_aircraft,"CD0",aircraft_idx)
    engine_type = InputValidate.get_value(df_aircraft,"Engine Type",aircraft_idx)
    number_of_engines = InputValidate.get_value(df_aircraft,"Number of Engines",aircraft_idx)

    # Initialise a variable to save TW data
    if engine_type == "Jet"
        TW = fill(0.0, length(WS), N_stages)
    else
        TW = fill(0.0 * u"m/s", length(WS), N_stages)
    end

    idx = 1

    # For each column of mission
    for col in (ncol(df_mission)-N_stages+1):ncol(df_mission)
        stage = df_mission[stage_row,col]

        if stage in ["Climb","Cruise","Loiter","Descend"]
            # Get the needed information
            h = uconvert(u"m", df_mission[h_row,col])
            a_max = uconvert(u"m/s^2", df_mission[a_max_row,col])
            engines_available = df_mission[engine_avail_row,col]
            gradient = df_mission[gradient_row,col]
            α = df_mission[α_row,col]
            β = df_mission[β_row,col]
            n = df_mission[n_row,col]

            if col in key_idx
                col_name = velocity_list[findfirst(==(col),velocity_list[:,"Saved Column"]),"Design Parameter"]
                V∞ = df_WS[!,col_name]
                V∞ = uconvert.(u"m/s", V∞)
            else
                V∞ = uconvert(u"m/s", df_mission[V∞_row,col])
            end
        
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

            TW[:,idx] = (α ./ β) .* (number_of_engines ./ engines_available) .* prop_constant .* (C_climb .+ C_speed .+ C_zerol .+ C_induc)
            name = DataFrames.names(df_mission)[col]
            df_WS[!,"TW $name"] = TW[:,idx]
        end
        idx += 1
    end

    return (TW, df_WS)
end

"""
    `second_seg_climb` - A function to obtain the second segment climb angle

    return the angles (OEI and AEO)
"""
function second_seg_climb(;df_aircraft::DataFrame,aircraft_idx::Int,V)
    V_V_AEO = InputValidate.get_value(df_aircraft,"Climb Rate AEO",aircraft_idx)
    V_V_OEI = InputValidate.get_value(df_aircraft,"Climb Rate OEI",aircraft_idx)
    γ_2_desc = InputValidate.get_value(df_aircraft,"Climb Gradient",aircraft_idx)
    N_engines = InputValidate.get_value(df_aircraft,"Number of Engines",aircraft_idx)

    γ_2_AEO = NaN
    γ_2_OEI = NaN

    if !ismissing(V_V_AEO)
        γ_2_AEO = asind(V_V_AEO / V)
    end

    if !ismissing(V_V_OEI)
        γ_2_OEI = asind(V_V_OEI / V)
    end

    # For FAR25
    if γ_2_desc == "FAR25"
        if N_engines == 2
            γ_2_OEI = atand(0.024)
        elseif N_engines == 3
            γ_2_OEI = atand(0.027)
        else
            γ_2_OEI = atand(0.030)
        end
    end

    return (γ_2_AEO, γ_2_OEI)
end

"""
    `takeoff_distance` - A function to obtain takeoff and initial climb TW or PW ratios

    return the TW or PW ratios
"""
function takeoff_distance(;df_WS,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame)
    # Find stages which are takeoff
    stage_idx = findfirst(==("Stage"),df_mission[:,1])
    col_idx = findall(==("Takeoff"), skipmissing(collect(df_mission[stage_idx, :])))

    WS = df_WS[:,1]
    
    ### Physical properties
    (ρ_0,_,_,_) = ISAdata(0*u"m") # Determine density
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant

    ### Aerodynamic coefficients
    # ASSUME C_L = 0 for takeoff, as it is effectively no lift (zero AoA)
    CL = fill(0.0, length(WS), 1)

    CL_max = InputValidate.get_value(df_aircraft,"Wing CLmax",aircraft_idx)
    CL_max_tof = InputValidate.get_value(df_aircraft,"Wing CLmax Takeoff Flaps",aircraft_idx)
    CL_max_total = CL_max + CL_max_tof

    CD0 = InputValidate.get_value(df_aircraft,"CD0",aircraft_idx)
    CD0_tof = InputValidate.get_value(df_aircraft,"Wing CD0 Takeoff Penalty",aircraft_idx)
    CD0_ldg = InputValidate.get_value(df_aircraft,"Wing CD0 Landing Gear Penalty",aircraft_idx)
    CD0_total = CD0 + CD0_tof + CD0_ldg

    AR = InputValidate.get_value(df_aircraft,"Wing AR",aircraft_idx)
    e = InputValidate.get_value(df_aircraft,"Oswald Efficiency",aircraft_idx)
    e_tof = InputValidate.get_value(df_aircraft,"Wing Oswald Efficiency Takeoff Penalty",aircraft_idx)
    e_ldg = InputValidate.get_value(df_aircraft,"Wing Oswald Efficiency Landing Gear Penalty",aircraft_idx)
    e_total = e + e_tof + e_ldg

    ### Runway conditions
    μ = InputValidate.get_value(df_aircraft,"Takeoff Friction",aircraft_idx)

    # If missing, means the takeoff friction has not been specified by the airworthiness requirements
    if ismissing(μ)
        # Take the runway friction
        df_friction = ParseData.get_data(data = "Runway")
    end

    ### Takeoff factors
    V_to_factor = InputValidate.get_value(df_aircraft,"Takeoff Velocity Relative to Stall",aircraft_idx)
    takeoff_distance_factor = InputValidate.get_value(df_aircraft,"Takeoff Distance Multiplier",aircraft_idx)
    engine_type = InputValidate.get_value(df_aircraft,"Engine Type",aircraft_idx)

    ### BFL Conditions -> Using Imperial units
    g_imp = uconvert(u"ft/s^2", g) # Gravitational acceleration constant
    obstacle_height_imp = uconvert(u"ft", InputValidate.get_value(df_aircraft,"Takeoff Obstacle",aircraft_idx))
    V_second_seg_factor = InputValidate.get_value(df_aircraft,"Climbing Velocity Relative to Stall",aircraft_idx)
    WS_imp = uconvert.(u"lbf/ft^2", WS)

    # Initialise Takeoff TW 
    if engine_type == "Jet"
        TW_ground_roll = fill(0.0, length(WS), length(col_idx))
        TW_climb_AEO = fill(0.0, length(WS), length(col_idx))
        TW_climb_OEI = fill(0.0, length(WS), length(col_idx))
        TW_BFL = fill(0.0, length(WS), length(col_idx))
    else
        TW_ground_roll = fill(0.0*u"m/s", length(WS), length(col_idx))
        TW_climb_AEO = fill(0.0*u"m/s", length(WS), length(col_idx))
        TW_climb_OEI = fill(0.0*u"m/s", length(WS), length(col_idx))
        TW_BFL = fill(0.0*u"m/s", length(WS), length(col_idx))
    end

    idx = 1

    for col in col_idx
        α = InputValidate.get_value(df_mission,"α",col)
        β = InputValidate.get_value(df_mission,"Engine Efficiency Scaling",col)
        ρ = InputValidate.get_value(df_mission,"ρ",col)
        ρ_imp = uconvert(u"slug/ft^3", ρ)
        σ = InputValidate.get_value(df_mission,"σ",col)
        takeoff_distance = InputValidate.get_value(df_mission,"Distance",col) * takeoff_distance_factor

        if ismissing(μ)
            runway = InputValidate.get_value(df_mission,"Runway",col)
            μ = InputValidate.get_value(df_friction,runway,"Friction")

            if ismissing(μ)
                throw(ArugmentError("Runway condition $runway cannot be found! Please check whether the runway condition is specified in the assumptions"))
            end
        end

        # Calculate ground roll TW ratio
        V_stall = sqrt.((2 .*α .* WS) ./ (ρ .* CL_max_total))
        V_LOF = V_stall .* V_to_factor
        K_A = ρ_0 ./ (2 * WS) .* (μ .* CL .- CD0_total .- ((CL_max_total .^ 2) ./ (pi * AR * e_total))) # Calculate K_A from ground roll equation
        TW_ground_roll_save = μ .+ ((K_A .* V_LOF .^2) ./ (exp.(2 .* g .* K_A .* takeoff_distance).-1))

        # Second Segment Climb Information
        takeoff_distance_imp = uconvert(u"ft", takeoff_distance)
        V_second_seg = V_stall .* V_second_seg_factor

        ### Determine whether the input is power or thrust 
        if engine_type == "Jet"
            λ_bpr = InputValidate.get_value(df_aircraft,"Engine Bypass Ratio",aircraft_idx)
            k_e = 0.75 * ((5 + λ_bpr) / (4 + λ_bpr))
        else
            prop_disc_load = uconvert(u"hp/ft^2", InputValidate.get_value(df_aircraft,"Engine Propeller Disc Load",aircraft_idx))
            prop_type = InputValidate.get_value(df_aircraft,"Engine Propeller Type",aircraft_idx)

            if prop_type == "Constant Speed"
                k_c = 1
            elseif prop_type == "Fixed Pitch"
                k_c = 0.8
            else
                throw("Type of Propeller must be either constant speed or fixed pitch!")
            end

            k_e = 5.75 * k_c * cbrt(σ / prop_disc_load)
        end

        # Get Second Segment Climb and Calculate TW ratios
        (γ_2_AEO, γ_2_OEI) = second_seg_climb(df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,V=V_second_seg)

        # BFL Calculations
        γ_2 = γ_2_OEI
        Δγ2 = γ_2 - γ_2_OEI # Assume worse case, so no 
        μ_prime = 0.01 * CL_max_total + μ
        CL_second_seg = CL_max_total ./ (V_second_seg_factor .^ 2)
        CD_second_seg = CD0_total .+ (1 ./ (pi*AR*e)).*CL_second_seg.^2 
        CL_CD_second_sg = CL_second_seg./CD_second_seg
        TW_climb_AEO_save = sind(γ_2_AEO) +  1 / CL_CD_second_sg
        TW_climb_OEI_save = sind(γ_2_OEI) + 1 / CL_CD_second_sg # ASSUME NO WINDMILLING
        TW_BFL_save = (α ./ (β * k_e)) .* ((1 ./ ((1.159 * (takeoff_distance_imp - 655u"ft" / sqrt(σ)) * (1 + 2.3 * Δγ2)) ./ (((WS_imp .* α) ./ (ρ_imp * g_imp * CL_second_seg) .+ obstacle_height_imp)) .- 2.7)) .+ μ_prime)

        # Convert to power to weight ratio for prop
        if engine_type in ["Propeller","Turboprop"]
            # Ground Roll TW
            μ_prop = InputValidate.get_value(df_aircraft,"Engine Propeller Efficiency",aircraft_idx)
            TW_ground_roll[:,idx] = TW_ground_roll_save .* V_LOF ./ μ_prop

            # Climb TW
            TW_climb_AEO[:,idx] .= TW_climb_AEO_save * u"m/s"
            TW_climb_OEI[:,idx] .= TW_climb_OEI_save * u"m/s"

            # BFL
            TW_BFL[:,idx] = ustrip.(upreferred.(TW_BFL_save)) * u"m/s"
        else
            TW_ground_roll[:,idx] = TW_ground_roll_save
            TW_climb_AEO[:,idx] .= TW_climb_AEO_save
            TW_climb_OEI[:,idx] .= TW_climb_OEI_save
            TW_BFL[:,idx] = TW_BFL_save
        end

        # Save into df
        name = DataFrames.names(df_mission)[col]
        df_WS[!,"TW Ground Roll $name"] = TW_ground_roll[:,idx]
        df_WS[!,"TW Climb AEO $name"] = TW_climb_AEO[:,idx]
        df_WS[!,"TW Climb OEI $name"] = TW_climb_OEI[:,idx]
        df_WS[!,"TW BFL $name"] = TW_BFL[:,idx]

        idx += 1
    end

    return (TW_ground_roll,TW_climb_AEO,TW_climb_OEI,TW_BFL,df_WS)
end

"""
    `get_max_WS` - A function to get the maximum wing loading from the landing conditions

    return the TW
"""
function get_max_WS(;df_aircraft::DataFrame, aircraft_idx::Int,df_mission::DataFrame)
    row_idx = findfirst(==("Stage"),df_mission[:,1])
    landing_col = findall(==("Landing"),skipmissing(collect(df_mission[row_idx, :])))

    
    K_l = 1/InputValidate.get_value(df_aircraft,"Landing Distance Multiplier",aircraft_idx) # Multiplier for landing distance (smaller)
    S_a = uconvert(u"m", InputValidate.get_value(df_aircraft,"Glide Slope",aircraft_idx)) # Glide Slope extra distance
    CL_max = InputValidate.get_value(df_aircraft,"Wing CLmax",aircraft_idx)
    ΔCL_max = InputValidate.get_value(df_aircraft,"Wing CLmax Landing Flaps",aircraft_idx)
    CL_max_land = CL_max + ΔCL_max
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant

    WS_max = 999999

    # For each landing column
    for col in landing_col
        # Get landing distance
        LDA = uconvert(u"m", InputValidate.get_value(df_mission,"Distance",col))
        α = InputValidate.get_value(df_mission,"Max_Mass_Ratio",col)
        reverser = InputValidate.get_value(df_mission,"Reverser",col)
        σ = InputValidate.get_value(df_mission,"σ",col)

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
function quick_constraint(;df_WS::DataFrame,velocity_list::DataFrame,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    # Get the power/thrust from the conditions except for takeoff
    (TW, _) = point_perf(df_WS=df_WS,velocity_list=velocity_list,df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages)

    # Takeoff TW
    (TW_ground_roll,TW_climb_AEO,TW_climb_OEI,TW_BFL) = takeoff_distance(df_WS=df_WS,df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission)

    # Concat the TW data
    TW_to_max = vcat(TW_ground_roll,TW_climb_AEO,TW_climb_OEI)

    BFL_consider = InputValidate.get_value(df_aircraft,"BFL Consideration",aircraft_idx)

    # If BFL has to be considered, then concatenate it
    if BFL_consider == true
        TW_to_max = vcat(TW_to_max,TW_BFL)
    end

    TW_to_max = TW_to_max[.!isnan.(TW_to_max)] # Remove all NaN

    # Get maxmimum values
    TW_to_max = maximum(TW_to_max)
    TW_other_max = maximum(TW)
    TW_max = max(TW_to_max,TW_other_max)

    # Return the results from the search
    return TW_max
end

"""
    `full_constraint` - A function to get all T/W

    return the saved WS df
"""
function full_constraint(;df_WS::DataFrame,velocity_list::DataFrame,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    # Get the power/thrust from the conditions except for takeoff
    (_, df_WS) = point_perf(df_WS=df_WS,velocity_list=velocity_list,df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages)

    # Takeoff TW
    (_,_,_,_,df_WS) = takeoff_distance(df_WS=df_WS,df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission)

    # Return the results from the search
    return df_WS
end

end