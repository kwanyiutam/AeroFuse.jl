# Main module for operation-centered design
module OperationOptimisation

# Initialise packages used
using ISAData
using Unitful
using CSV
using DataFrames
using QuasiMonteCarlo

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the ParseData module
include("parsedata.jl")

# Include the InputValidate module
include("inputvalidate.jl")

# include the CostModel module
include("costmodel.jl")

# include the InitialSizing module
include("initialsizing.jl")

# include the ConstraintDiagram module
include("constraintdiagram.jl")

# include the AircraftOptimisation module
include("aircraftoptimisation.jl")

"""
    `design_init` - A function which initiates the design variables

    Returns the relevant dataframes required for design
"""
function design_init(;mission_path :: String = "Default", aircraft_path :: String = "Default", payload_path :: String = "Default", cost_path :: String = "Default")
    # Initiate a Design Parameter DataFrame
    design_param = DataFrame(
        "Design Parameter" => String[],
        "Saved Location" => String[],
        "Saved Row" => Int64[],
        "Saved Column" => Int64[],
        "Type" => String[],
        "Continuous" => Bool[],
        "Lower Bound" => Float64[],
        "Upper Bound" => Float64[]
    )

    ### These contains the aircrafts to compare
    if aircraft_path == "Default"
        df_aircraft = ParseData.get_data(data = "SampleAircraft");
    else
        # Try and find the dataset with the matching "data" ID 
        df_aircraft = try
            CSV.read(aircraft_path, DataFrame)
        catch e
            throw(ArgumentError("Cannot find the relevant dataset in $aircraft_path !"))
        end
    end

    (df_aircraft, design_param, N_aircraft) = InputValidate.validate_data(df = df_aircraft,design_param = design_param, data_type = "Aircraft")

    ### These contains the payload parameters
    if payload_path == "Default"
        df_payload = ParseData.get_data(data = "SamplePayload");
    else
        df_payload = try
            CSV.read(payload_path, DataFrame)
        catch e
            throw(ArgumentError("Cannot find the relevant dataset in $payload_path !"))
        end
    end

    (df_payload, design_param, N_payload) = InputValidate.validate_data(df = df_payload,design_param = design_param, data_type = "Payload")
    
    ### These contains the baseline mission to compare
    if mission_path == "Default"
        df_mission = ParseData.get_data(data = "SampleMission");
    else
        df_mission = try
            CSV.read(mission_path, DataFrame)
        catch e
            throw(ArgumentError("Cannot find the relevant dataset in $mission_path !"))
        end
    end

    (df_mission, design_param, N_stages) = InputValidate.validate_data(df = df_mission,design_param = design_param, data_type = "Mission")

    df_cost = CostModel.cost_init()

    return (design_param,df_aircraft,N_aircraft,df_mission,N_stages,df_payload,N_payload,df_cost)
end

"""
    `sampling_points` - A function to sample points in a design space

"""
function sampling_points(;N_sample::Int, design_param::DataFrame, design_type::String="Design")
    design_row_idx = findall(==(design_type), design_param[:,"Type"])

    if isempty(design_row_idx)
        return DataFrame()
    end

    design_names = design_param[design_row_idx,"Design Parameter"]
    design_lb = design_param[design_row_idx,"Lower Bound"]
    design_ub = design_param[design_row_idx,"Upper Bound"]
    continuous_check = design_param[design_row_idx,"Continuous"]
    
	sample_points = QuasiMonteCarlo.sample(N_sample, design_lb, design_ub, LatinHypercubeSample())

    df_samples = DataFrame(sample_points',design_names)

    for idx in 1:ncol(df_samples)
        if continuous_check[idx] == false
            df_samples[:,idx] = round.(Int, df_samples[:,idx])
        end
    end

    return df_samples
end

"""
    `WS_init` - A function which initialise the WS dataframe to store key information

"""
function WS_init(;WS,velocity_list::DataFrame,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame)
    # Create a Wing loading dataframe, filled with just WS for now
    df_WS = DataFrame(WS = WS)
    (ρ_0,_,_,_) = ISAdata(0*u"m")

    # Fill it with the design velocity as two other columns
    vel_unit = AeroUnits.convert_to_unit(df_mission[findfirst(==("Velocity"),df_mission[:,1]),"Units"])
    for i in 1:nrow(velocity_list)
        df_WS[!, velocity_list[i,"Design Parameter"]] = fill(0.0*vel_unit, length(WS))
    end

    df_WS[!, "MTOW"]  = fill(0.0*u"kg", length(WS))
    df_WS[!, "Empty Weight"]  = fill(0.0*u"kg", length(WS))
    df_WS[!, "Fuel Weight"]  = fill(0.0*u"kg", length(WS))
    df_WS[!, "Total Cost"]  = fill(0.0, length(WS))
    df_WS[!, "Fuel Cost"]  = fill(0.0, length(WS))
    df_WS[!, "Cabin Crew Cost"]  = fill(0.0, length(WS))
    df_WS[!, "Flight Crew Cost"]  = fill(0.0, length(WS))

    # Define key variables
    AR = df_aircraft[findfirst(==("Wing AR"),df_aircraft[:,1]),aircraft_idx]
    e = df_aircraft[findfirst(==("Oswald Efficiency"),df_aircraft[:,1]),aircraft_idx]
    CD0 = df_aircraft[findfirst(==("CD0"),df_aircraft[:,1]),aircraft_idx]
    CL_max = df_aircraft[findfirst(==("Wing CLmax"),df_aircraft[:,1]),aircraft_idx]
    
    # Update the V_imd column at MTOW
    df_WS[!, "V_imd_MTOW"] = sqrt.(2*WS./ρ_0).*((1/(pi*AR*e*CD0)).^0.25)
    df_WS[!, "V_stall_MTOW_clean"] = sqrt.(2*WS./(ρ_0*CL_max))

    return df_WS
end

function update_design(;WS_max, TW_max, df_aircraft::DataFrame, N_aircraft::Int, aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    MTOW = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx)
    AR = InputValidate.get_value(df_aircraft,"Wing AR",aircraft_idx)

    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant

    if MTOW > 150.0*1000*u"kg"
        cost_model = "Full_A"
    elseif MTOW > 51.0*1000*u"kg"
        cost_model = "Full_B"
    elseif MTOW > 15.0*1000*u"kg"
        cost_model = "Full_C"
    else
        cost_model = "Full_D"
    end

    # Calculate the values
    Sref = MTOW * g / WS_max
    Tmax = MTOW * g * TW_max
    bref = sqrt(AR * Sref)

    # Calculate the maximum velocity and mach number
    vel_row = findfirst(==("Velocity"),df_mission[:,1])
    alt_row = findfirst(==("Altitude"),df_mission[:,1])
    df_mission_filtered = df_mission[:,ncol(df_mission)-N_stages+1:ncol(df_mission)]
    max_idx = argmax(skipmissing(collect(df_mission_filtered[vel_row,:])))
    V_max = df_mission_filtered[vel_row,max_idx]
    (rho, P, _, _) = ISAdata(uconvert(u"m",df_mission_filtered[alt_row,max_idx]))
    c_sound = upreferred(sqrt(1.4 * P / rho))
    M_max = uconvert(u"m/s",V_max) / c_sound

    # Calculate the design cruise range
    dist_row = findfirst(==("Distance"),df_mission[:,1])
    range_max = maximum(skipmissing(collect(df_mission_filtered[dist_row,:])))

    # Calculate maximum landing weight
    mlw_row = findfirst(==("Max_Mass_Ratio"),df_mission[:,1])
    α_land_idx = argmax(skipmissing(collect(df_mission_filtered[mlw_row,:])))
    α_land = df_mission_filtered[mlw_row,α_land_idx]
    MLW = α_land * MTOW

    # Get landing stall speed
    ρ_row = findfirst(==("ρ"),df_mission[:,1])
    ρ = df_mission_filtered[ρ_row,α_land_idx]
    CL_max_total = InputValidate.get_value(df_aircraft, "Wing CLmax", aircraft_idx) + InputValidate.get_value(df_aircraft, "Wing CLmax Landing Flaps", aircraft_idx)
    V_s_land = sqrt.((2 .* α_land .* WS_max) ./ (ρ .* CL_max_total))

    # Calculate block time
    block_time = sum(df_mission[findfirst(==("Duration"),df_mission[:,1]),ncol(df_mission)-N_stages+1:ncol(df_mission)])

    # Append to the aircraft data!
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="WS_max",value=WS_max,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Wing Area",value=Sref,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Wing Span",value=bref,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Tmax",value=Tmax,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Cost Model",value=cost_model,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Block Time",value=block_time,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Operating Velocity",value=V_max,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Operating Mach Number",value=M_max,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Maximum Range",value=range_max,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="MLW",value=MLW,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Landing Stall Speed",value=V_s_land,N_config=N_aircraft,col=aircraft_idx)

    return df_aircraft
end

"""
    `mission_design_flow` - A function which runs the design workflow
"""
function mission_design_flow(;design_param::DataFrame,velocity_list::DataFrame,df_pert::DataFrame,df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,N_payload::Int,df_cost::DataFrame,df_save::DataFrame)
    for payload_idx in (ncol(df_payload)-N_payload+1):ncol(df_payload)
        (df_aircraft,df_mission,df_payload) = InitialSizing.update_init_design(df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,payload_idx=payload_idx) 
            
        # Assume WS_max is restricted by landing!
        WS_max = ConstraintDiagram.get_max_WS(df_aircraft = df_aircraft, aircraft_idx = aircraft_idx,df_mission = df_mission)

        # Only need to optimise at WS_max, because that is where the minimum cost is likely at (Quick mode)
        df_WS = WS_init(WS=[WS_max],velocity_list=velocity_list,df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission)
        
        # Get the minimum cost velocity and their respective design
        (df_WS, df_mission) = InitialSizing.velocity_optimise_main(df_WS = df_WS,velocity_list = velocity_list,df_aircraft=df_aircraft,N_aircraft=N_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,payload_idx=payload_idx,df_cost=df_cost)

        # Get the T/W ratio required
        TW_max = ConstraintDiagram.quick_constraint(WS_max=WS_max,df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages)

        # Initialise the aircraft design
        df_aircraft = update_design(WS_max = WS_max, TW_max = TW_max, df_aircraft = df_aircraft, N_aircraft = N_aircraft, aircraft_idx = aircraft_idx,df_mission=df_mission,N_stages=N_stages)

        # Run the aircraft optimisation / perturbations
        df_save = AircraftOptimisation.aircraft_optimisation_start(design_param = design_param,df_pert = df_pert,df_aircraft=df_aircraft,N_aircraft=N_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,payload_idx=payload_idx,df_cost=df_cost,df_save=df_save)
    end
    
    return df_save
end


"""
    `design_start` - A function which runs the design workflow

"""
function design_start(;design_param::DataFrame,df_design::DataFrame,df_pert::DataFrame, df_aircraft::DataFrame,N_aircraft::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,N_payload::Int,df_cost::DataFrame,data_save::Vector{String})
    # Get the list of design-specific parameters
    design_list = design_param[findall(==("Design"),design_param[:,:Type]),:]

    # Get a list of velocity
    velocity_list = design_param[findall(value -> occursin(Regex("(?i)Velocity"),value),design_param[:,:Type]),:]

    # Check that the df_design inputs are complete: It should have four columns for each design-specific parameter
    if sort(lowercase.(DataFrames.names(df_design))) != sort(lowercase.(design_list[:,"Design Parameter"]))
        throw(ArgumentError("Invalid df_design input, some design paramters did not match / not complete."))
    end
    
    # Define key variables needed for initial sizing
    df_aircraft = InitialSizing.define_aircraft_properties(df_aircraft=df_aircraft,N_aircraft=N_aircraft)

    # Default entries
    data_save_def = ["Aircraft Index", "Payload Index"]
    data_save_def = vcat(data_save_def,design_param[:,"Design Parameter"],"Total Cost","Cost Breakdown","Wing Mesh","HT Mesh","VT Mesh","Fuselage Shape","Engine Shape")
    data_save_def = vcat(data_save_def, data_save)
    df_save = DataFrame(NamedTuple{Tuple(Symbol.(data_save_def))}(Tuple(Any[] for _ in data_save_def)))

    # For each aircraft
    for aircraft_idx in (ncol(df_aircraft)-N_aircraft+1):ncol(df_aircraft)
        if nrow(df_design) == 0
            # Directly run the mission design flow
            df_save = mission_design_flow(design_param=design_param,df_pert=df_pert,df_aircraft=df_aircraft,N_aircraft=N_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,df_cost=df_cost,df_save=df_save)
        else
            # For each specified design point
            for row in 1:nrow(df_design)
                ### Update properties based on the design parameters
                # For each design parameter
                for param in 1:ncol(df_design)
                    param_name = DataFrames.names(df_design)[param]
                    param_value = df_design[row,param]
            
                    # Update the datasets to include the new parameters
                    (df_aircraft,df_mission,df_payload) = InputValidate.update_df_with_design(design_list=design_list,parameter=param_name,value=param_value,df_aircraft=df_aircraft,df_mission=df_mission,df_payload=df_payload)
                end

                # Run the mission design flow
                df_save = mission_design_flow(design_param=design_param,velocity_list=velocity_list,df_pert=df_pert,df_aircraft=df_aircraft,N_aircraft=N_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,df_cost=df_cost,df_save=df_save)
            end
        end
    end

    return df_save
end

end
