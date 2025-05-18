module AircraftOptimisation

# Initialise packages used
using DataFrames
using Unitful

# Include the InputValidate module
include("inputvalidate.jl")

# Include the AircraftAero module
include("aircraftaero.jl")

# include the CostModel module
include("costmodel.jl")

"""
    `fuselage_sizing` - A function which sizes the fuselage to the appropriate size

    Returns the updated dataframes
"""
function fuselage_sizing(;df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_payload::DataFrame,N_payload::Int,payload_col::Int)
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
   
    # Save diameter into aircraft
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Diameter",value=D_eff,N_config=N_aircraft,col=aircraft_idx)

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

    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Cabin Length",value=cabin_length,N_config=N_aircraft,col=aircraft_idx)

    # Calculate ratios and length properties for the aircraft
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Cabin Fineness Ratio",value=cabin_length/D_eff,N_config=N_aircraft,col=aircraft_idx)
    nose_fineness = df_aircraft[findfirst(==("Nose Fineness Ratio"),df_aircraft[:,1]),aircraft_idx]
    nose_length = nose_fineness*D_eff
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Nose Length",value=nose_length,N_config=N_aircraft,col=aircraft_idx)
    tail_fineness = df_aircraft[findfirst(==("Afterbody Fineness Ratio"),df_aircraft[:,1]),aircraft_idx]
    tail_length = tail_fineness*D_eff
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Afterbody Length",value=tail_length,N_config=N_aircraft,col=aircraft_idx)
    fuselage_length = nose_length + tail_length + cabin_length
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Fuselage Length",value=fuselage_length,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Fuselage Fineness Ratio",value=fuselage_length/D_eff,N_config=N_aircraft,col=aircraft_idx)
    
    return (df_aircraft, df_payload)
end

"""
    `powerplant_sizing` - A function which sizes the powerplant's dimensions and weight (Rubber Engine)
"""
function powerplant_sizing(;df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int)
    engine_type = InputValidate.get_value(df_aircraft,"Engine Type",aircraft_idx)

    if engine_type == "Jet"
        bpr = InputValidate.get_value(df_aircraft,"Engine Bypass Ratio",aircraft_idx)
        Tmax = ustrip(InputValidate.get_value(df_aircraft,"Tmax",aircraft_idx)) / 1000.0 # In kN
        M = ustrip(InputValidate.get_value(df_aircraft,"Operating Mach Number",aircraft_idx))
        
        # SFC for high bypass ratio
        if bpr > 6
            @warn "Bypass Ratio too big for this regression model! SFC will assume to be same as initial sizing, while weight and dimensions will be assumed using BPR = 6"
            SFC_cruise = InputValidate.get_value(df_aircraft,"SFC_cruise",aircraft_idx)
            SFC_loiter = InputValidate.get_value(df_aircraft,"SFC_loiter",aircraft_idx)
            bpr = 6
        else
            SFC_cruise = 25*exp(-0.05*bpr) / 1000 *u"g / N / s"
            SFC_loiter = SFC_cruise * (InputValidate.get_value(df_aircraft,"SFC_loiter",aircraft_idx) / InputValidate.get_value(df_aircraft,"SFC_cruise",aircraft_idx))
        end

        # Afterburner or non-afterburner engines
        if bpr < 1
            weight = 11.1*Tmax^1.1*M^0.25*exp(-0.81*bpr)*u"kg"
            length = 0.68*Tmax^0.4*M^0.2*u"m"
            diameter = 0.11*Tmax^0.5*exp(0.04*bpr)*u"m"
            SFC_cruise = 30*exp(-0.186*bpr) / 1000 *u"g / N / s"
            SFC_loiter = SFC_cruise * (InputValidate.get_value(df_aircraft,"SFC_loiter",aircraft_idx) / InputValidate.get_value(df_aircraft,"SFC_cruise",aircraft_idx))
        else
            weight = 14.7*Tmax^1.1*exp(-0.045*bpr)*u"kg"
            length = 0.49*Tmax^0.4*M^0.2*u"m"
            diameter = 0.15*Tmax^0.5*exp(0.04*bpr)*u"m"
        end
    else
        # Get the power
        Pmax = ustrip(InputValidate.get_value(df_aircraft,"Tmax",aircraft_idx)) / 1000.0 # In kW

        # No data on SFC, just use initial sizing data
        SFC_cruise = InputValidate.get_value(df_aircraft,"SFC_cruise",aircraft_idx)
        SFC_loiter = InputValidate.get_value(df_aircraft,"SFC_loiter",aircraft_idx)

        if engine_type == "Turboprop"
            if Pmax < 300
                @warn "Power is below the applicable range, the data is now being extrapolated. Proceed with caution"
            elseif Pmax > 3728
                @warn "Power is above the applicable range, the data is now being extrapolated. Proceed with caution"
            end

            weight = 0.96*Pmax^0.803*u"kg"
            length = 0.12*Pmax^0.373*u"m"
            diameter = 0.25*Pmax^0.120*u"m"
        else
            weight = 3.12*Pmax^0.780*u"kg"
            length = 0.11*Pmax^0.424*u"m"
            diameter = 0.8*u"m"
        end
    end

    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Engine Weight Estimate",value=weight,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Engine Length Estimate",value=length,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Engine Diameter Estimate",value=diameter,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="SFC Cruise Estimate",value=SFC_cruise,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="SFC Loiter Estimate",value=SFC_loiter,N_config=N_aircraft,col=aircraft_idx)

    return df_aircraft
end

"""
    `aircraft_design_flow` - A function which runs the design workflow
"""
function aircraft_design_flow(;opt_list::DataFrame,df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,N_payload::Int,payload_idx::Int,df_cost_model::DataFrame)
    # Update fuselage information based on the number of passengers
    (df_aircraft, df_payload) = fuselage_sizing(df_aircraft = df_aircraft,N_aircraft = N_aircraft, aircraft_idx = aircraft_idx, df_payload=df_payload, N_payload=N_payload, payload_col=payload_idx)

    # More refined update for weights
    df_aircraft = powerplant_sizing(df_aircraft = df_aircraft,N_aircraft=N_aircraft,aircraft_idx=aircraft_idx)

    # Given the wing, fuselage information, calcualate aerodynamic properties + generate mesh
    (df_aircraft, df_mission) = AircraftAero.run_aero_analysis(df_aircraft=df_aircraft,N_aircraft=N_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages)
end

"""
    `aircraft_optimisation_start` - A function which starts the optimisation process for design (with perturbation first)
"""
function aircraft_optimisation_start(;design_param::DataFrame,df_pert::DataFrame,df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,N_payload::Int,payload_idx::Int,df_cost::DataFrame)
    # List of perturbation / optimisation parameters
    pert_list = design_param[findall(==("Perturbations"),design_param[:,:Type]),:]
    opt_list = design_param[findall(==("Optimise"),design_param[:,:Type]),:]

    # Check that the df_pert inputs are complete: It should have four columns for each design-specific parameter
    if sort(lowercase.(DataFrames.names(df_pert))) != sort(lowercase.(pert_list[:,"Design Parameter"]))
        throw(ArgumentError("Invalid df_pert input, some design paramters did not match / not complete."))
    end

    cost_model = df_aircraft[findfirst(==("Cost Model"),df_aircraft[:,1]),aircraft_idx]
    df_cost_model = CostModel.cost_model(df_cost = df_cost, model = cost_model)

    if nrow(pert_list) == 0
        # Directly run the mission design flow
        aircraft_design_flow(opt_list=opt_list,df_aircraft=df_aircraft,N_aircraft=N_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,payload_idx=payload_idx,df_cost_model=df_cost_model)
    else
        # For each specified design point
        for row in 1:nrow(df_pert)
            ### Update properties based on the design parameters
            # For each design parameter
            for param in 1:ncol(df_pert)
                param_name = DataFrames.names(df_pert)[param]
                param_value = df_pert[row,param]
        
                # Update the datasets to include the new parameters
                (df_aircraft,df_mission,df_payload) = InputValidate.update_df_with_design(design_list=pert_list,parameter=param_name,value=param_value,df_aircraft=df_aircraft,df_mission=df_mission,df_payload=df_payload)
            end

            # Run the aircraft design evaluation flow
            aircraft_design_flow(opt_list=opt_list,df_aircraft=df_aircraft,N_aircraft=N_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,payload_idx=payload_idx,df_cost_model=df_cost_model)
        end
    end
end

end