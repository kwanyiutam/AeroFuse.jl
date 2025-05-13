module InitialSizing
# Initialise packages used
using CSV
using DataFrames
using LsqFit
using Unitful
using ISAData

using Optimization
using ForwardDiff

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the ParseData module
include("parsedata.jl")

# Include the InputValidate module
include("inputvalidate.jl")

# include the CostModel module
include("costmodel.jl")

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
            units = u"g / N / s"
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
        W_cargo = df_payload[findfirst(==("Cargo"),df_payload[:,1]),col]

        # Calculate total payload weight
        W_payload = W_ppl + W_bag + W_cargo + W_drop
        df_payload = InputValidate.df_update_or_append(df=df_payload,label="Payload Weight",value=W_payload,N_config=N_payload,col=col)
    end

    (ρ_0,_,_,_) = ISAdata(0*u"m")

    for col in (ncol(df_mission)-N_stages+1):ncol(df_mission)
        (ρ,_,_,_) = ISAdata(df_mission[findfirst(==("Altitude"),df_mission[:,1]),col])
        σ = ρ / ρ_0

        df_mission = InputValidate.df_update_or_append(df=df_mission,label="σ",value=σ,N_config=N_stages,col=col)
    end

    return (df_aircraft,df_mission,df_payload)
end

function weight_converge(;W_0::Float64, W_payload::Float64, fuel_ratio::Float64, A::Float64, C::Float64)
    W_0 = ustrip(uconvert(u"kg", W_0))
    W_payload = ustrip(uconvert(u"kg", W_payload))

    W_e = NaN
    W_0_prev = 0.0
    max_iter = 50
    iter = 1

    while abs(W_0-W_0_prev) > 0.01
        # Empty weight Estimate
        W_e = A*(W_0^C)
        W_0_prev = W_0

        # New guess for MTOW
        W_0 = (W_payload) / (1 - fuel_ratio - W_e)

        iter = iter + 1

        # Checks to ensure the loop will not run forever
        if W_0 < 0
            W_0 = NaN
            W_e = NaN
            @warn "MTOW returned a negative value, so this is invalid!"
            break
        elseif iter >= max_iter
            @warn "Maximum iteration $max_iter was reached! Proceed with caution as the results may not be fully converged."
            break
        end
    end

    return(W_0, W_e)
end

function fuel_weight_fractions(;stage::String,engine_type::String,V_md::Float64,V_bar::Float64, LD_max::Float64,SFC_cruise::Float64,SFC_loiter::Float64, endurance::Float64 = 0.0, range::Float64 = 0.0)
    if engine_type != "Jet"
        SFC = SFC * V_bar * V_md
    end

    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant


    if stage == "Takeoff"
        fraction = 0.97
    elseif stage == "Climb"
        fraction = 0.985
    elseif stage == "Cruise"
        @assert range > 0 "Range must be given!"
        # Convert units
        range = uconvert(u"m", range)
        fraction = exp(-(endurance * SFC_cruise * g * (V_bar + V_bar^-3)) / (LD_max * V_md * 2))
    elseif stage == "Loiter"
        @assert endurance > 0 "Endurance must be given!"
        # Convert units
        endurance = uconvert(u"s", endurance)
        fraction = exp(-(endurance * SFC_loiter * g * (V_bar^2 + V_bar^-2)) / (LD_max * 2))
    else
        # assume mostly idle thrust, so descend condition
        fraction = 0.995
    end

    return fraction
end

function fuel_ratio(;df_WS::DataFrame,WS_idx::Int,min_fuel_vel::DataFrame,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    # Trapped fuel ratio (2%)
    trapped_fuel = 1.02

    # Get the minimum fuel columns
    min_fuel_stages = min_fuel_vel[:,"Saved Column"]

    # Append new row/update row so that takeoff is mass ratio of missing (to be added in the for-loop)
    df_mission = InputValidate.df_update_or_append(df=df_mission,label="Mass_Ratio",value=missing,N_config=N_stages,col=ncol(df_mission)-N_stages+1)
    # Alpha is defined as the ending weight fractions
    df_mission = InputValidate.df_update_or_append(df=df_mission,label="Cumulative_Mass_Ratio",value=missing,N_config=N_stages,col=ncol(df_mission)-N_stages+1)

    V_imd_MTOW = df_WS[WS_idx, "V_imd_MTOW"]

    engine_type = df_aircraft[findfirst(==("Engine Type"),df_aircraft[:,1]),aircraft_idx]
    LD_max = df_aircraft[findfirst(==("LD_max"),df_aircraft[:,1]),aircraft_idx]
    SFC_cruise = df_aircraft[findfirst(==("SFC_cruise"),df_aircraft[:,1]),aircraft_idx]
    SFC_loiter = df_aircraft[findfirst(==("SFC_loiter"),df_aircraft[:,1]),aircraft_idx]

    # For each mission stages
    for col in (ncol(df_mission)-N_stages+1):ncol(df_mission)
        stage = df_mission[findfirst(==("Stage"),df_mission[:,1]),col]
        σ = df_mission[findfirst(==("σ"),df_mission[:,1]),col]
        
        # Define alpha - the cumulative weight fractions - for calculations
        if col == ncol(df_mission)-N_stages+1
            α = 1
        else
            α = df_mission[findfirst(==("Cumulative_Mass_Ratio"),df_mission[:,1]),col-1]
        end

        # Get the minimum drag speed
        V_md = V_imd_MTOW*sqrt(α/σ)
        print(V_md)

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

        endurnace = 0
        range = 0

        if stage == "Cruise"
            range = df_mission[findfirst(==("Distance"),df_mission[:,1]),col]
        elseif stage == "Loiter"
            endurance = df_mission[findfirst(==("Duration"),df_mission[:,1]),col]
        end

        fraction = fuel_weight_fractions(stage=stage,engine_type=engine_type,V_md=V_md,V_bar=V_bar,LD_max=LD_max,SFC_cruise=SFC_cruise,SFC_loiter=SFC_loiter,endurance=endurance,range=range)
        df_mission = InputValidate.df_update_or_append(df=df_mission,label="Mass_Ratio",value=fraction,N_config=N_stages,col=col)

        if col == ncol(df_mission)-N_stages+1
            α = fraction
        else
            α = df_mission[findfirst(==("Cumulative_Mass_Ratio"),df_mission[:,1]),col-1] * fraction
        end

        df_mission = InputValidate.df_update_or_append(df=df_mission,label="Cumulative_Mass_Ratio",value=α,N_config=N_stages,col=col)
    end

    # Final fuel ratio * trapped fuel
    final_fuel_ratio = (1-df_mission[findfirst(==("Cumulative_Mass_Ratio"),df_mission[:,1]),ncol(df_mission)])*trapped_fuel

    return (final_fuel_ratio,df_mission)
end

function opt_basic_cost(x, p)
    # Expand the parameters of p
    df_WS = p.df_WS
    WS_idx = p.WS_idx
    aircraft_idx = p.aircraft_idx
    N_stages = p.N_stages
    payload_idx = p.payload_idx
    min_cost_vel = p.min_cost_vel
    min_fuel_vel = p.min_fuel_vel
    df_aircraft = p.df_aircraft
    df_mission = p.df_mission
    df_payload = p.df_payload
    df_cost = p.df_cost

    # Load guesses of x into the dataframe
    (df_aircraft,df_mission,df_payload) = InputValidate.update_df_with_design(design_list=min_cost_vel,parameter=min_cost_vel[:,"Design Parameter"],value=x,df_aircraft=df_aircraft,df_mission=df_mission,df_payload=df_payload)

    # Calculate fuel ratios
    (final_fuel_ratio,df_mission) = fuel_ratio(df_WS=df_WS,WS_idx=WS_idx,min_fuel_vel=min_fuel_vel,df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages)

    # Get key values
    W_payload = df_payload[findfirst(==("Payload Weight"),df_payload[:,1]),payload_idx]
    A = df_aircraft[findfirst(==("A"),df_aircraft[:,1]),aircraft_idx]
    C = df_aircraft[findfirst(==("C"),v[:,1]),aircraft_idx]

    # Converge weight, starting with a random weight
    (MTOW, empty_weight) = weight_converge(W_0 = 100000*u"kg", W_payload = W_payload, fuel_ratio=final_fuel_ratio, A=A, C=C)

    (total_cost, _) = cost_calc(df_cost = df_cost, model = "Basic")
end

function velocity_optimise_main(;df_WS::DataFrame,velocity_list::DataFrame,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,payload_idx::Int,df_cost::DataFrame)
    #V_imd = sqrt.(2*α.*WS./ρ_0).*((1/(pi*aircraft.AR*aircraft.e*aircraft.CD0)).^0.25)
    #V_md = V_imd ./ sqrt(σ)

    # Get the minimum drag speed for each wing loading
    #df_WS[:,"V_md"] = df_WS[:,"WS"]

    # Get the velocities
    min_fuel_vel = velocity_list[findall(value -> occursin(Regex("(?i)Minimum Fuel Velocity"),value),velocity_list[:,:Type]),:]
    min_cost_vel = velocity_list[findall(value -> occursin(Regex("(?i)Minimum Cost Velocity"),value),velocity_list[:,:Type]),:]

    # Optimisation Parameters
    x0 = [250.0, 250.0] # Initial guesses --> Just random numbers for now
    lb = min_cost_vel[:,"Lower Bound"]
    ub = min_cost_vel[:,"Upper Bound"]

    df_cost_model = CostModel.cost_model(df_cost = df_cost, model = "Basic")
     
    for idx in 1:nrow(df_WS)
        # Hyperparameters for optimisation
        p = (df_WS = df_WS, WS_idx = idx, min_cost_vel = min_cost_vel, min_fuel_vel = min_fuel_vel, df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,payload_idx=payload_idx)

        # Problem construction
        #optprob = OptimizationFunction(opt_basic_cost, Optimization.AutoForwardDiff())
        #prob = OptimizationProblem(optprob, x0, p, lb = lb, ub = ub)

        ## Choose optimizer
        #solution = solve(prob, BFGS())

        #print(solution)
    end

    return (df_WS, df_mission)
end

end