module InitialSizing
# Initialise packages used
using CSV
using DataFrames
using LsqFit
using Unitful
using ISAData

#using Optimization
#using OptimizationBBO
#using ForwardDiff

#using NonlinearSolve

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the ParseData module
include("parsedata.jl")

# Include the InputValidate module
include("inputvalidate.jl")

# include the CostModel module
include("costmodel.jl")

function fuselage_sizing(;df_payload::DataFrame,N_payload::Int,payload_col::Int)
    # Get the volume required for cargo bay (baggage + cargo + droppable)
    W_baggage = df_payload[findfirst(==("Baggage Total Weight"),df_payload[:,1]),payload_col]
    W_cargo = df_payload[findfirst(==("Cargo Weight"),df_payload[:,1]),payload_col]
    density_baggage = df_payload[findfirst(==("Baggage Packing Density"),df_payload[:,1]),payload_col]
    density_cargo = df_payload[findfirst(==("Cargo Packing Density"),df_payload[:,1]),payload_col]
    V_baggage = uconvert(u"m^3", W_baggage/density_baggage)
    V_cargo = uconvert(u"m^3", W_cargo/density_cargo)
    V_storage = V_baggage + V_cargo

    # Fuselage Wall thickness (one side, so multiply by two for diameter)
    wall_thick = uconvert(u"m", df_payload[findfirst(==("Fuselage Wall Thickness"),df_payload[:,1]),payload_col])

    # Get the number of passengers and initiate different number of rows
    N_passenger = df_payload[findfirst(==("Passengers"),df_payload[:,1]),payload_col]
    seat_pitch = uconvert(u"m", df_payload[findfirst(==("Seat Pitch"),df_payload[:,1]),payload_col])
    seat_abreast = df_payload[findfirst(==("Seat Abreast"),df_payload[:,1]),payload_col]
    N_rows = ceil(N_passenger / seat_abreast)
    seat_length = seat_pitch * N_rows

    # Properties for height
    floor_thick = df_payload[findfirst(==("Floor Thickness Relative to Diameter"),df_payload[:,1]),payload_col]
    cargo_height = uconvert(u"m", df_payload[findfirst(==("Cargo Bay Height"),df_payload[:,1]),payload_col])
    head_room = uconvert(u"m", df_payload[findfirst(==("Head Room"),df_payload[:,1]),payload_col])
    extra_height = uconvert(u"m", df_payload[findfirst(==("Extra Height"),df_payload[:,1]),payload_col])
    approx_height = (head_room+cargo_height+(2*wall_thick)+extra_height)*(1+floor_thick) # Just an approximation, likely need further refinement

    # Properties for width
    aisle_width = uconvert(u"m", df_payload[findfirst(==("Aisle Width"),df_payload[:,1]),payload_col])
    cargo_width = uconvert(u"m", df_payload[findfirst(==("Cargo Bay Width"),df_payload[:,1]),payload_col])
    extra_width = uconvert(u"m", df_payload[findfirst(==("Extra Width"),df_payload[:,1]),payload_col])

    # Get seat widths of different configs
    seat_width_1 = uconvert(u"m", df_payload[findfirst(==("Seat Width 1 Pax"),df_payload[:,1]),payload_col])
    seat_width_2 = uconvert(u"m", df_payload[findfirst(==("Seat Width 2 Pax"),df_payload[:,1]),payload_col])
    seat_width_3 = uconvert(u"m", df_payload[findfirst(==("Seat Width 3 Pax"),df_payload[:,1]),payload_col])

    # Start filling from three across then down
    N_seat_3 = floor(seat_abreast / 3)
    N_seat_2 = floor((seat_abreast-N_seat_3*3) / 2)
    N_seat_1 = seat_abreast-N_seat_3*3-N_seat_2*2

    seat_width_total = seat_width_1*N_seat_1 + seat_width_2*N_seat_2 + seat_width_3*N_seat_3

    approx_width = max((seat_width_total+aisle_width),cargo_width) + 2*wall_thick + extra_width

    # Effective diameter
    D_eff = sqrt(approx_width*approx_height)

    df_payload = InputValidate.df_update_or_append(df=df_payload,label="Seat Rows",value=N_rows,N_config=N_payload,col=payload_col)
    df_payload = InputValidate.df_update_or_append(df=df_payload,label="Three Seats",value=N_seat_3,N_config=N_payload,col=payload_col)
    df_payload = InputValidate.df_update_or_append(df=df_payload,label="Two Seats",value=N_seat_2,N_config=N_payload,col=payload_col)
    df_payload = InputValidate.df_update_or_append(df=df_payload,label="One Seat",value=N_seat_1,N_config=N_payload,col=payload_col)
    df_payload = InputValidate.df_update_or_append(df=df_payload,label="Diameter",value=D_eff,N_config=N_payload,col=payload_col)

    # Length calculations
    cross_aisle_length = uconvert(u"m", df_payload[findfirst(==("Cross Aisle per 20"),df_payload[:,1]),payload_col])
    lavatories_length = uconvert(u"m", df_payload[findfirst(==("Lavatories per 50"),df_payload[:,1]),payload_col])
    galley_vol = uconvert(u"m^3", df_payload[findfirst(==("Galley per pax"),df_payload[:,1]),payload_col])

    # Get the actual length increase
    cross_aisle_length = cross_aisle_length * ceil(N_passenger / 20)
    lavatories_total_length = lavatories_length * ceil(N_passenger / 50)

    # Lavatories are square, and can be packed tighter (right now assume all lengthwise)
    N_lavatory_rows = ceil(D_eff / lavatories_total_length)
    lavatory_actual_length = N_lavatory_rows * lavatories_length

    # Approximate the galley length needed by assuming head room x seat width is all used for galley, just an approximation
    galley_length = (galley_vol * N_passenger) / (head_room * seat_width_total)

    # Get cabin length
    extra_length = uconvert(u"m", df_payload[findfirst(==("Extra Length"),df_payload[:,1]),payload_col])
    cabin_length = seat_length + cross_aisle_length + lavatory_actual_length + galley_length + extra_length
    
    # Check how much cargo space left is needed
    V_cargo_remain = V_storage - (cargo_height*cargo_width*cabin_length)

    # If still need to carry more cargo, assume a separate cargo space at the back of the cabin (extending the cabin length)
    if V_cargo_remain > 0.0 * u"m^3"
        # Assume the back cargo bay is effectively the same as the belly one
        cabin_length = V_cargo_remain / (cargo_width * cargo_height)
    end

    df_payload = InputValidate.df_update_or_append(df=df_payload,label="Cabin Length",value=cabin_length,N_config=N_payload,col=payload_col)

    return df_payload
end

"""
    `define_aircraft_properties` - A function which defines aircraft properties needed for initial sizing

    Returns the updated dataframes
"""
function define_aircraft_properties(;df_aircraft::DataFrame,N_aircraft::Int)
    # For each aircraft configuration column
    for col in (ncol(df_aircraft)-N_aircraft+1):ncol(df_aircraft)
        # Calculate LD_max
        KLD = df_aircraft[findfirst(==("KLD"),df_aircraft[:,1]),col]
        AR = df_aircraft[findfirst(==("Wing AR"),df_aircraft[:,1]),col]
        SwetSref = df_aircraft[findfirst(==("Swet/Sref"),df_aircraft[:,1]),col]
        e = df_aircraft[findfirst(==("Oswald Efficiency"),df_aircraft[:,1]),col]

        LD_max = KLD*sqrt(AR/SwetSref)
        df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="LD_max",value=LD_max,N_config=N_aircraft,col=col)

        # Calculate CD0
        CD0 = pi*AR*e / ((2*LD_max)^2)
        df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="CD0",value=CD0,N_config=N_aircraft,col=col)

        engine_type = df_aircraft[findfirst(==("Engine Type"),df_aircraft[:,1]),col]
        # Update SFC
        if engine_type == "Jet"
            units = u"g / N / s"
            λ_bpr = df_aircraft[findfirst(==("Engine Bypass Ratio"),df_aircraft[:,1]),col]

            # High bypass ratio definition (Arbitrary for now)
            if λ_bpr > 5
                SFC_cruise = 14.1
                SFC_loiter = 11.3
            # Low bypass ratio definition
            elseif λ_bpr > 0
                SFC_cruise = 22.7
                SFC_loiter = 19.8
            # Pure turbojet
            else
                SFC_cruise = 25.5
                SFC_loiter = 22.7
            end

            SFC_cruise = SFC_cruise / 1000 * units # Convert from mg to g
            SFC_loiter = SFC_loiter / 1000 * units # Convert from mg to g
        else
            units = u"g / W / s"
            prop_type = df_aircraft[findfirst(==("Engine Propeller Type"),df_aircraft[:,1]),col]
            μ_prop = df_aircraft[findfirst(==("Engine Propeller Efficiency"),df_aircraft[:,1]),col]

            # Turboprop
            if engine_type == "Turboprop"
                SFC_cruise = 0.085
                SFC_loiter = 0.101
            # Low bypass ratio definition
            elseif prop_type == "Fixed Pitch"
                SFC_cruise = 0.068
                SFC_loiter = 0.085
            # Pure turbojet
            elseif prop_type == "Constant Speed"
                SFC_cruise = 0.068
                SFC_loiter = 0.085
            else
                throw(ErrorArgument("Propeller Type $prop_type is not valid!"))
            end

            # IF PROP, WOULD NEED TO MULTIPLY BY VELOCITY!
            SFC_cruise = SFC_cruise / 1000 * units / μ_prop # Convert from mg to g
            SFC_loiter = SFC_loiter / 1000 * units / μ_prop # Convert from mg to g
        end

        df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="SFC_cruise",value=SFC_cruise,N_config=N_aircraft,col=col)
        df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="SFC_loiter",value=SFC_loiter,N_config=N_aircraft,col=col)
    end
    
    return df_aircraft
end

"""
    `update_init_design` - A function which updates properties needed for initial sizing

    Returns the updated dataframes
"""
function update_init_design(;df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,N_payload::Int)
    W_drop = sum(df_mission[findfirst(==("Payload_Drop"),df_mission[:,1]),(ncol(df_mission)-N_stages+1):ncol(df_mission)])
    
    # For each payload configuration column
    for col in (ncol(df_payload)-N_payload+1):ncol(df_payload)
        # Get number of passenger
        N_pax = df_payload[findfirst(==("Passengers"),df_payload[:,1]),col]

        # Get number of flight crews, using minimum flight crew for now, and append to the payload row
        N_fcrew = df_aircraft[findfirst(==("Min Flight Crew"),df_aircraft[:,1]),aircraft_idx]
        df_payload = InputValidate.df_update_or_append(df=df_payload,label="Flight Crew",value=N_fcrew,N_config=N_payload,col=col)

        # Get number of cabin crews, and append to the payload row
        N_ccrew = 0

        # 1 cabin crew every 50 passenger approximation
        if N_pax > 10
            N_ccrew = ceil(Int64, N_pax / 50)
        end

        df_payload = InputValidate.df_update_or_append(df=df_payload,label="Cabin Crew",value=N_ccrew,N_config=N_payload,col=col)

        # Get passenger + crew weight
        W_per_ppl = df_payload[findfirst(==("Human Weight"),df_payload[:,1]),col]
        W_ppl = (N_pax+N_fcrew+N_ccrew)*W_per_ppl

        # Calculate Baggage weight per person
        W_bag_per_ppl = df_payload[findfirst(==("Baggage Weight"),df_payload[:,1]),col]
        W_bag = (N_pax+N_fcrew+N_ccrew)*W_bag_per_ppl

        # Calculate Cargo weight
        W_cargo = df_payload[findfirst(==("Cargo Weight"),df_payload[:,1]),col]

        # Calculate total payload weight
        W_payload = W_ppl + W_bag + W_cargo + W_drop
        df_payload = InputValidate.df_update_or_append(df=df_payload,label="Baggage Total Weight",value=W_bag,N_config=N_payload,col=col)
        df_payload = InputValidate.df_update_or_append(df=df_payload,label="Droppable Total Weight",value=W_bag,N_config=N_payload,col=col)
        df_payload = InputValidate.df_update_or_append(df=df_payload,label="Payload Weight",value=W_payload,N_config=N_payload,col=col)

        # Update fuselage information based on the number of passenger
        df_payload = fuselage_sizing(df_payload=df_payload, N_payload=N_payload, payload_col=col)
    end

    (ρ_0,_,_,_) = ISAdata(0*u"m")

    for col in (ncol(df_mission)-N_stages+1):ncol(df_mission)
        (ρ,_,_,_) = ISAdata(df_mission[findfirst(==("Altitude"),df_mission[:,1]),col])
        σ = ρ / ρ_0

        df_mission = InputValidate.df_update_or_append(df=df_mission,label="ρ",value=ρ,N_config=N_stages,col=col)
        df_mission = InputValidate.df_update_or_append(df=df_mission,label="σ",value=σ,N_config=N_stages,col=col)
    end

    return (df_aircraft,df_mission,df_payload)
end

"""
    `weight_converge` - A function to iterate for the MTOW and empty weight for initial sizing

    Returns the MTOW and empty weights
"""
function weight_converge(;W_0_init, W_payload, fuel_ratio, A, C)
    W_0_init = ustrip(uconvert(u"kg", W_0_init))
    W_payload = ustrip(uconvert(u"kg", W_payload))

    W_e_W_0 = NaN
    W_0_prev = 0.0
    W_0 = W_0_init
    max_iter = 100
    iter = 1

    while abs(W_0-W_0_prev) > 0.01
        # Empty weight Estimate
        W_e_W_0 = A*(W_0^C)
        W_0_prev = W_0

        # New guess for MTOW
        W_0 = (W_payload) / (1 - fuel_ratio - W_e_W_0)

        iter = iter + 1

        # Checks to ensure the loop will not run forever
        if W_0 < 0
            # Return NaN
            return (NaN, NaN)
            #@warn "MTOW returned a negative value, so this is invalid!"
            break
        elseif iter >= max_iter
            #@warn "Maximum iteration $max_iter was reached! Proceed with caution as the results may not be fully converged."
            break
        end
    end

    W_0 = W_0*u"kg"
    W_e = W_e_W_0 * W_0

    return(W_0, W_e)
end

"""
    `fuel_weight_fractions` - A function to obtain the fuel fractions at different stages of flight

    Returns the fuel fractions
"""
function fuel_weight_fractions(;stage::String,engine_type::String,V_md,V_bar,LD_max,SFC_cruise,SFC_loiter,endurance = 0.0, range= 0.0)
    if engine_type != "Jet"
        SFC_cruise = upreferred(SFC_cruise * V_bar * V_md)
        SFC_loiter = upreferred(SFC_loiter * V_bar * V_md)
    end

    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant

    if stage == "Takeoff"
        fraction = 0.97
    elseif stage == "Climb"
        fraction = 0.985
    elseif stage == "Cruise"
        @assert ustrip(range) > 0.0 "Range must be given!"
        # Convert units
        range = uconvert(u"m", range)
        fraction = exp(-(range * SFC_cruise * g * (V_bar + V_bar^-3)) / (LD_max * V_md * 2))
    elseif stage == "Loiter"
        @assert ustrip(endurance) > 0.0 "Endurance must be given!"
        # Convert units
        endurance = uconvert(u"s", endurance)
        fraction = exp(-(endurance * SFC_loiter * g * (V_bar^2 + V_bar^-2)) / (LD_max * 2))
    else
        # assume mostly idle thrust, so descend condition
        fraction = 0.995
    end

    return fraction
end

"""
    `fuel_ratio` - A function to obtain the final fuel fraction after considering all stages

    Returns the fuel fractions and dataframe
"""
function fuel_ratio(;df_WS::DataFrame,WS_idx::Int,min_fuel_vel::DataFrame,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,save_info::Bool = false)
    # Trapped fuel ratio (2%)
    trapped_fuel = 1.02

    # Get the minimum fuel columns
    min_fuel_stages = min_fuel_vel[:,"Saved Column"]
    V_imd_MTOW = df_WS[WS_idx, "V_imd_MTOW"]

    engine_type = df_aircraft[findfirst(==("Engine Type"),df_aircraft[:,1]),aircraft_idx]
    LD_max = df_aircraft[findfirst(==("LD_max"),df_aircraft[:,1]),aircraft_idx]
    SFC_cruise = df_aircraft[findfirst(==("SFC_cruise"),df_aircraft[:,1]),aircraft_idx]
    SFC_loiter = df_aircraft[findfirst(==("SFC_loiter"),df_aircraft[:,1]),aircraft_idx]
    cumulative_weight_fraction = 1

    # For each mission stages
    for col in (ncol(df_mission)-N_stages+1):ncol(df_mission)
        stage = df_mission[findfirst(==("Stage"),df_mission[:,1]),col]
        σ = df_mission[findfirst(==("σ"),df_mission[:,1]),col]

        # Get the minimum drag speed
        V_md = V_imd_MTOW*sqrt(cumulative_weight_fraction/σ)

        if col in min_fuel_stages
            if (engine_type == "Jet" && stage == "Cruise") || (engine_type in ["Turboprop","Propeller"] && stage == "Loiter")
                V = V_md * 3^(1/4)
            else
                V = V_md
            end

            df_mission = InputValidate.df_update_or_append(df=df_mission,label="Velocity",value=V,N_config=N_stages,col=col)
        else
            V = df_mission[findfirst(==("Velocity"),df_mission[:,1]),col]
        end

        V_bar = V / V_md

        endurance = 0.0
        range = 0.0

        if stage == "Cruise"
            range = df_mission[findfirst(==("Distance"),df_mission[:,1]),col]
        elseif stage == "Loiter"
            endurance = df_mission[findfirst(==("Duration"),df_mission[:,1]),col]
        end

        fraction = fuel_weight_fractions(stage=stage,engine_type=String(engine_type),V_md=V_md,V_bar=V_bar,LD_max=LD_max,SFC_cruise=SFC_cruise,SFC_loiter=SFC_loiter,endurance=endurance,range=range)

        if col == ncol(df_mission)-N_stages+1
            cumulative_weight_fraction = fraction
        else
            cumulative_weight_fraction *= fraction
        end

        if save_info == true
            df_mission = InputValidate.df_update_or_append(df=df_mission,label="Fuel Fraction",value=fraction,N_config=N_stages,col=col)
        end

        if stage != "Loiter"
            duration = 15*u"minute" # For Takeoff and Landing

            if stage in ["Climb", "Descend"]
                duration = 30*u"minute"
            elseif stage == "Cruise"
                # Convert to time
                duration = uconvert(u"minute", range / V)
            end

            df_mission = InputValidate.df_update_or_append(df=df_mission,label="Duration",value=duration,N_config=N_stages,col=col)
        end
    end

    # Final fuel ratio * trapped fuel
    final_fuel_ratio = (1-cumulative_weight_fraction)*trapped_fuel

    return (final_fuel_ratio,df_mission)
end

"""
    `basic_cost_calc` - A function to calculate the basic operational costs

    Returns the costs and relevant dataframes
"""
function basic_cost_calc(x, p, save_info::Bool = false)
    # Expand the parameters of p
    df_WS = p.df_WS
    WS_idx = p.WS_idx
    min_cost_vel = p.min_cost_vel
    min_fuel_vel = p.min_fuel_vel
    df_aircraft = p.df_aircraft
    N_aircraft = p.N_aircraft
    aircraft_idx = p.aircraft_idx
    df_mission = p.df_mission
    N_stages = p.N_stages
    df_payload = p.df_payload
    payload_idx = p.payload_idx
    df_cost_model = p.df_cost_model

    # Load guesses of x into the dataframe
    for i in eachindex(x)
        (df_aircraft,df_mission,df_payload) = InputValidate.update_df_with_design(design_list=min_cost_vel,parameter=min_cost_vel[i,"Design Parameter"],value=x[i],df_aircraft=df_aircraft,df_mission=df_mission,df_payload=df_payload)
    end

    # Calculate fuel ratios
    (final_fuel_ratio,df_mission) = fuel_ratio(df_WS=df_WS,WS_idx=WS_idx,min_fuel_vel=min_fuel_vel,df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,save_info=save_info)

    # Get key values
    W_payload = df_payload[findfirst(==("Payload Weight"),df_payload[:,1]),payload_idx]
    A = df_aircraft[findfirst(==("A"),df_aircraft[:,1]),aircraft_idx]
    C = df_aircraft[findfirst(==("C"),df_aircraft[:,1]),aircraft_idx]

    # Converge weight, starting with a random weight
    (MTOW, empty_weight) = weight_converge(W_0_init = 1000000*u"kg", W_payload = W_payload, fuel_ratio=final_fuel_ratio, A=A, C=C)

    if isnan(MTOW)
        return (NaN, NaN, df_aircraft, df_mission, df_payload)
    end

    # Calculate the fuel weight
    fuel_weight = final_fuel_ratio * MTOW

    # Append the weights into the aircraft data
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="MTOW",value=MTOW,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Empty Weight",value=empty_weight,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Fuel Weight",value=fuel_weight,N_config=N_aircraft,col=aircraft_idx)

    (total_cost, df_cost_breakdown) = CostModel.cost_calc(df_cost_model = df_cost_model, df_aircraft = df_aircraft, aircraft_idx = aircraft_idx, df_mission = df_mission, N_stages = N_stages, df_payload = df_payload, payload_idx = payload_idx)

    return (total_cost, df_cost_breakdown, df_aircraft, df_mission, df_payload)
end

"""
    `opt_basic_cost` - A function to calculate the basic operational costs

    Only returns the total cost (relevant for optimisations if needed)
"""
function opt_basic_cost(x, p)
    (total_cost, _, _, _, _) = basic_cost_calc(x, p)

    return total_cost
end

"""
    `velocity_optimise_main` - A main function to optimise velocity (for min cost)

    Rerturn updated values for df_WS and df_mission
"""
function velocity_optimise_main(;df_WS::DataFrame,velocity_list::DataFrame,df_aircraft::DataFrame,N_aircraft::Int, aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,payload_idx::Int,df_cost::DataFrame)
    # Get the velocities
    min_fuel_vel = velocity_list[findall(value -> occursin(Regex("(?i)Minimum Fuel Velocity"),value),velocity_list[:,:Type]),:]
    min_cost_vel = velocity_list[findall(value -> occursin(Regex("(?i)Minimum Cost Velocity"),value),velocity_list[:,:Type]),:]

    # Optimisation Parameters
    #x0 = fill(300.0, nrow(min_cost_vel)) # Initial guesses --> Just random numbers for now
    lb = min_cost_vel[:,"Lower Bound"]
    ub = min_cost_vel[:,"Upper Bound"]

    df_cost_model = CostModel.cost_model(df_cost = df_cost, model = "Basic")
     
    for idx in 1:nrow(df_WS)
        # Hyperparameters for optimisation
        p = (df_WS = df_WS, WS_idx = idx, min_cost_vel = min_cost_vel, min_fuel_vel = min_fuel_vel, df_aircraft=df_aircraft, N_aircraft = N_aircraft, aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,payload_idx=payload_idx,df_cost_model=df_cost_model)

        ### Old Optimisation attempt, did not go well!
        # Problem construction
        #optprob = OptimizationFunction(opt_basic_cost)
        #prob = OptimizationProblem(optprob, x0, p, lb = lb, ub = ub)

        ## Choose optimizer
        #solution = solve(prob, OptimizationBBO.BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 100000, maxtime = 1000.0)

        #if idx == nrow(df_WS)
            #print(solution)
        #end

        # Cannot be smaller than stall speed
        V_stall = df_WS[idx, "V_stall_MTOW_clean"]

        V_save = fill(300.0, nrow(min_cost_vel))

        for V_idx in 1:nrow(min_cost_vel)
            # Get the minimum and maximum velocities
            V_min = min(ustrip(V_stall),lb[V_idx])
            V_max = ub[V_idx]

            V_iterate = LinRange(V_min, V_max, 500)

            min_V_idx = 1
            min_cost = maxintfloat(Float64)

            for idx in eachindex(V_iterate)
                V_save[V_idx] = V_iterate[idx]

                (total_cost, _, _, _, _) = basic_cost_calc(V_save, p)

                if isnan(total_cost)
                    continue
                elseif min_cost > total_cost
                    min_cost = total_cost
                    min_V_idx = idx
                end
            end

            # Save the minimum speed and continue to iterate
            global V_save[V_idx] = V_iterate[min_V_idx]
        end

        # Final iteration with the minimum cost speed optimised
        (total_cost, df_cost_breakdown, df_aircraft, df_mission, df_payload) = basic_cost_calc(V_save, p, true)

        for vel in 1:nrow(velocity_list)
            name = velocity_list[vel, "Design Parameter"]
            row = velocity_list[vel, "Saved Row"]
            col = velocity_list[vel, "Saved Column"]
            df_WS[idx,name] = df_mission[row,col]
        end

        MTOW = df_aircraft[findfirst(==("MTOW"),df_aircraft[:,1]),aircraft_idx]

        df_WS[idx, "MTOW"]  = MTOW
        df_WS[idx, "Empty Weight"]  = df_aircraft[findfirst(==("Empty Weight"),df_aircraft[:,1]),aircraft_idx]
        df_WS[idx, "Fuel Weight"]  = df_aircraft[findfirst(==("Fuel Weight"),df_aircraft[:,1]),aircraft_idx]
        df_WS[idx, "Total Cost"]  = total_cost
        df_WS[idx, "Fuel Cost"]  = df_cost_breakdown[findfirst(==("Fuel"),df_cost_breakdown[:,1]),"Cost (USD)"]
        df_WS[idx, "Cabin Crew Cost"]  = df_cost_breakdown[findfirst(==("Cabin Crew"),df_cost_breakdown[:,1]),"Cost (USD)"]
        df_WS[idx, "Flight Crew Cost"]  = df_cost_breakdown[findfirst(==("Flight Crew"),df_cost_breakdown[:,1]),"Cost (USD)"]

        α = 1
        for col in ncol(df_mission)-N_stages+1:ncol(df_mission)
            fuel_fraction = df_mission[findfirst(==("Fuel Fraction"),df_mission[:,1]),col]
            dropped_payload = df_mission[findfirst(==("Payload_Drop"),df_mission[:,1]),col]
            drop_ratio = 1 - (dropped_payload / MTOW)
            α = α*fuel_fraction*drop_ratio
            df_mission = InputValidate.df_update_or_append(df=df_mission,label="α",value=α,N_config=N_stages,col=col)
        end
    end

    return (df_WS, df_mission)
end

end