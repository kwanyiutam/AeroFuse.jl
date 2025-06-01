module MultiMission

# Initialise packages used
using CSV
using DataFrames
using Unitful
using AeroFuse
using Roots
using Statistics
using Plots
using Interpolations
using LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font, guidefont = font(14,plot_font), tickfont = font(12,plot_font), legendfont = font(8,plot_font), titlefont = font(16,plot_font), colorbar_titlefont = font(12,plot_font))

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the InputValidate module
include("inputvalidate.jl")

# Include the AircraftAero module
include("aircraftaero.jl")

# Include the CostModel module
include("costmodel.jl")

# include the InitialSizing module
include("initialsizing.jl")

"""
    `AFD_calcs` - A function to calculate average flight durations (AFD) for utilisation
"""
function AFD_calcs(df_ops::DataFrame,df_aircraft::DataFrame,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame)
    if nrow(df_ops) != 0
        # Get df_ops columns
        ops_col_names = DataFrames.names(df_ops)

        # Initiate key cruise columns
        key_cruise_col = []
        new_duration = []

        for col_name in ops_col_names
            ac_data = findfirst(==(col_name),df_aircraft[:,1])
            payload_data = findfirst(==(col_name),df_payload[:,1])
            mission_data = findfirst(==(col_name),DataFrames.names(df_mission))

            if !isnothing(ac_data)
                throw(ArgumentError("Sorry, $col_name is an aircraft property! It cannot be changed between missions!"))
            elseif !isnothing(payload_data)
                # Ignore payload data
                continue
            elseif !isnothing(mission_data)
                value = df_mission[:,mission_data]
                stage_idx = findfirst(==("Stage"),df_mission[:,1])
                stage = value[stage_idx]

                # If the mission change is related to cruise
                if stage == "Cruise"
                    # Get current velocity
                    V_idx = findfirst(==("Velocity"),df_mission[:,1])
                    V = value[V_idx]

                    # Get all current range
                    range = df_ops[:,col_name]

                    # Get original range
                    range_idx = findfirst(==("Distance"),df_mission[:,1])
                    old_range = value[range_idx]
                    range = range .* unit(old_range)

                    # Calculate block time using that velocity
                    t_B = uconvert.(u"m",range) ./ uconvert(u"m/s",V)
                    t_B = uconvert.(u"minute",t_B)

                    # The duration for this column would change
                    push!(key_cruise_col,mission_data)
                    push!(new_duration,mean(t_B))
                else
                    throw(ArgumentError("Sorry, $stage missions are not currently supported for operaitonal analysis!"))
                end
            else
                throw(ErrorException("Cannot find the data $col_name on any of the datasets"))
            end
        end

        key_idx = ncol(df_mission)-N_stages+1:ncol(df_mission)
        key_idx = setdiff(key_idx, key_cruise_col)
        AFD = sum(df_mission[findfirst(==("Duration"),df_mission[:,1]),key_idx]) + sum(new_duration)
    else
        AFD = sum(df_mission[findfirst(==("Duration"),df_mission[:,1]),ncol(df_mission)-N_stages+1:ncol(df_mission)])
    end

    return AFD
end

"""
    `multi_mission_analysis` - A function that runs multi_mission analysis
"""
function multi_mission_analysis(df_ops::DataFrame,df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_payload::DataFrame,N_payload::Int,payload_idx::Int,df_mission::DataFrame,N_stages::Int,df_cost_model::DataFrame)
    if nrow(df_ops) != 0
        # Do deep copies so it does not affect outer loop
        df_aircraft = deepcopy(df_aircraft)
        df_payload = deepcopy(df_payload)
        df_mission = deepcopy(df_mission)

        # Obtain key aircraft details
        wing_mesh = InputValidate.get_value(df_aircraft,"Wing Mesh",aircraft_idx)
        Sw = InputValidate.get_value(df_aircraft,"Wing Area",aircraft_idx)
        HT_mesh = InputValidate.get_value(df_aircraft,"HT Mesh",aircraft_idx)
        VT_mesh = InputValidate.get_value(df_aircraft,"VT Mesh",aircraft_idx)
        fuse = InputValidate.get_value(df_aircraft,"Fuselage Shape",aircraft_idx)
        eng_save = InputValidate.get_value(df_aircraft,"Engine Shape",aircraft_idx)
        CD_upsweep = InputValidate.get_value(df_aircraft,"CD Upsweep",aircraft_idx)

        # Trapped fuel ratio (2%)
        trapped_fuel = InputValidate.get_value(df_aircraft,"Trapped Fuel Ratio",aircraft_idx)

        # Get MTOW
        MTOW = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx)
        g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
        engine_type = InputValidate.get_value(df_aircraft,"Engine Type",aircraft_idx)

        # Get df_ops columns
        ops_col_names = DataFrames.names(df_ops)

        # Parses the data and see where to retrive / update data later
        ops_param = DataFrame(
            "Design Parameter" => String[],
            "Saved Location" => String[],
            "Saved Row" => Int64[],
            "Saved Column" => Int64[],
        )

        key_cruise_monitor = []

        # Get the maximum values from design
        max_fuel_weight = InputValidate.get_value(df_aircraft,"Fuel Weight",aircraft_idx)
        max_payload_weight = InputValidate.get_value(df_payload,"Payload Weight",payload_idx)
        empty_weight = InputValidate.get_value(df_aircraft,"Empty Weight",aircraft_idx)

        # Prepare for calculating off-design performances
        for col_name in ops_col_names
            ac_data = findfirst(==(col_name),df_aircraft[:,1])
            payload_data = findfirst(==(col_name),df_payload[:,1])
            mission_data = findfirst(==(col_name),DataFrames.names(df_mission))

            if !isnothing(ac_data)
                throw(ArgumentError("Sorry, $col_name is an aircraft property! It cannot be changed between missions!"))
            elseif !isnothing(payload_data)
                save_design = [col_name,"Payload",payload_data,payload_idx]
                push!(ops_param,save_design)
            elseif !isnothing(mission_data)
                value = df_mission[:,mission_data]
                stage_idx = findfirst(==("Stage"),df_mission[:,1])
                stage = value[stage_idx]

                # If the mission change is related to cruise
                if stage == "Cruise"
                    range_row = findfirst(==("Distance"),df_mission[:,1])
                    save_design = [col_name,"Mission",range_row,mission_data]

                    push!(ops_param,save_design)
                    append!(key_cruise_monitor,mission_data)

                    vel_idx = findfirst(==("Velocity"),df_mission[:,1])
                    ρ_idx = findfirst(==("ρ"),df_mission[:,1])
                    V = uconvert(u"m/s",value[vel_idx])
                    ρ = uconvert(u"kg/m^3",value[ρ_idx])

                    refs = References(
                        speed     = ustrip(V),
                        area      = ustrip(Sw),
                        density   = ustrip(ρ),
                        span      = span(wing_mesh),
                        chord     = mean_aerodynamic_chord(wing_mesh),
                        location  = mean_aerodynamic_center(wing_mesh)
                    )

                    # Get the aircraft's CD0
                    α0 = try
                        find_zero((-15,15)) do α
                            sys = AircraftAero.make_case(α, wing_mesh, HT_mesh, VT_mesh, refs)
                            0.0 - AircraftAero.get_forces(sys, wing_mesh, HT_mesh, VT_mesh, fuse, eng_save, CD_upsweep).CL
                        end
                    catch e
                        @error "ERROR: " exception=(e, catch_backtrace())

                        @warn "Case Failed! This row will be returning NaN..."
                        return (NaN, DataFrame(), DataFrame(), [])
                        #df_mission = InputValidate.df_update_or_append(df=df_mission,label="Vmd_MTOW",value=NaN*u"m/s",N_config=N_stages,col=mission_data)
                        #df_mission = InputValidate.df_update_or_append(df=df_mission,label="LDmax",value=NaN,N_config=N_stages,col=mission_data)

                        #continue
                    end

                    sys = AircraftAero.make_case(α0, wing_mesh, HT_mesh, VT_mesh, refs)
                    init = AircraftAero.get_forces(sys, wing_mesh, HT_mesh, VT_mesh, fuse, eng_save, CD_upsweep)
                    CD0 = init.CD
                    
                    #CD_cruise = InputValidate.get_value(df_aircraft,"CD_cruise",aircraft_idx)
                    #CL_cruise = InputValidate.get_value(df_aircraft,"CL_cruise",aircraft_idx)
                    e = InputValidate.get_value(df_aircraft,"Oswald Efficiency",aircraft_idx)
                    AR = InputValidate.get_value(df_aircraft,"Wing AR",aircraft_idx)

                    #e_estimate = CL_cruise^2 / ((CD_cruise - CD0)*pi*AR)

                    Vmd = sqrt(2*MTOW*g / (ρ*Sw))*(1/(pi*AR*e*CD0))^0.25

                    LDmax = 0.5*sqrt((pi*AR*e)/CD0)
    
                    df_mission = InputValidate.df_update_or_append(df=df_mission,label="Vmd_MTOW",value=Vmd,N_config=N_stages,col=mission_data)
                    df_mission = InputValidate.df_update_or_append(df=df_mission,label="LDmax",value=LDmax,N_config=N_stages,col=mission_data)
                else
                    throw(ArgumentError("Sorry, $stage missions are not currently supported!"))
                end
            else
                throw(ErrorException("Cannot find the data $col_name on any of the datasets"))
            end
        end

        # Obtain Payload-Range Diagram Details
        idx_save = 1

        save_payload_range = []
        for col in key_cruise_monitor
            # Calculate different points on the payload-range diagram
            # 1 - Max Payload, MTOW
            # 2 - Max Fuel, MTOW
            # 3 - Max Fuel, max range
            fuel_weight_carried = [MTOW - max_payload_weight - empty_weight, max_fuel_weight, max_fuel_weight]
            payload_weight_carried = [max_payload_weight, MTOW - max_fuel_weight - empty_weight, 0.0*u"kg"]

            takeoff_weight = fuel_weight_carried .+ payload_weight_carried .+ empty_weight

            Vmd_mtow = uconvert(u"m/s",InputValidate.get_value(df_mission,"Vmd_MTOW",col))
            LD_max = InputValidate.get_value(df_mission,"LDmax",col)

            # Get the fuel fractions of the other columns except for the current one
            key_cols = setdiff(ncol(df_mission)-N_stages+1:ncol(df_mission), col)

            other_frac = prod(df_mission[findfirst(==("Fuel Fraction"),df_mission[:,1]),key_cols])

            # Get the new fuel fraction
            new_fuel_frac = (1 .- (fuel_weight_carried ./ takeoff_weight ./ trapped_fuel)) ./ other_frac

            # Get the previous weight fraction
            α = InputValidate.get_value(df_mission,"α",col-1)
            V_md = Vmd_mtow.*sqrt.(α)

            SFC_loiter = InputValidate.get_value(df_aircraft,"SFC_loiter",aircraft_idx)
            SFC_cruise = InputValidate.get_value(df_aircraft,"SFC_cruise",aircraft_idx)

            # Cruise velocity, take min fuel to get the payload-range diagram
            # SFC references
            if engine_type != "Jet"
                V_bar = 1.0 / (3.0^(0.25))
                V = uconvert.(u"m/s", V_md * V_bar)
                SFC_loiter = upreferred(SFC_loiter .* V)
                SFC_cruise = upreferred(SFC_cruise .* V)
            else
                V_bar = 3.0^(0.25)
                V = uconvert.(u"m/s", V_md * V_bar)
            end

            R = upreferred.(-(LD_max .* V_md .* 2 .* log.(new_fuel_frac))./ (SFC_cruise .* g .* (V_bar .+ V_bar.^-3)))

            # Add one more stage at the start:
            # When range = 0, can still carry max payload
            R = uconvert.(u"km",vcat(0.0*u"km",R))
            payload_weight_carried = uconvert.(u"kg",vcat(max_payload_weight,payload_weight_carried))

            duration = R ./ V
            push!(save_payload_range, [R , payload_weight_carried, fill(V, length(R)), duration])

            idx_save += 1
        end

        # Info to be saved
        total_cost_save = fill(NaN, nrow(df_ops))
        df_breakdown_save = DataFrame()
        V_optimal_save = fill(NaN*u"knots", nrow(df_ops))
        range_save = fill(NaN*u"km", nrow(df_ops))
        payload_save = fill(NaN*u"kg", nrow(df_ops))
        fuel_weight_save = fill(NaN*u"kg", nrow(df_ops))

        # Initiate an interpolation for payload cap (NOTE: since max fuel weight leads to MTOW, the 2nd and 3rd index is the same, so it is overwritten)
        idx = [1,2,4]

        intp = Interpolations.linear_interpolation(save_payload_range[1][1][idx], save_payload_range[1][2][idx], extrapolation_bc = Flat())
        intp_V = Interpolations.linear_interpolation(save_payload_range[1][1][idx], save_payload_range[1][3][idx], extrapolation_bc = Flat())
        intp_duration = Interpolations.linear_interpolation(save_payload_range[1][1][idx], save_payload_range[1][4][idx], extrapolation_bc = Flat())

        # For each row of operations
        for ops_idx in 1:nrow(df_ops)
            # Update each data
            for param in 1:ncol(df_ops)
                param_name = DataFrames.names(df_ops)[param]
                param_value = df_ops[ops_idx,param]
        
                # Update the datasets to include the new parameters
                (df_aircraft,df_mission,df_payload) = InputValidate.update_df_with_design(design_list=ops_param,parameter=param_name,value=param_value,df_aircraft=df_aircraft,df_mission=df_mission,df_payload=df_payload)
            end

            # Update the design (with new passengers number etc.)
            (df_aircraft,df_mission,df_payload) = InitialSizing.update_init_design(df_aircraft=df_aircraft,aircraft_idx=aircraft_idx,df_mission=df_mission,N_stages=N_stages,df_payload=df_payload,N_payload=N_payload,payload_idx=payload_idx)

            payload_save[ops_idx] = uconvert(u"kg", InputValidate.get_value(df_payload,"Payload Weight",payload_idx))

            # Initialise a range of velocities for optimisation
            V_iter = LinRange(100*u"knots", 500*u"knots", 3)
            #TODO: This currently does not support multiple cruise missions, because that would require V_iter to be different

            # Estimation for weight remaining
            α = fill(1.0, length(V_iter))

            # Save updated durations
            duration_save = fill(0.0*u"minute", length(V_iter), length(key_cruise_monitor))
            duration_save_old = fill(0.0*u"minute", 1, length(key_cruise_monitor))
            save_idx = 1

            for col in ncol(df_mission)-N_stages+1:ncol(df_mission)
                if col in key_cruise_monitor
                    Vmd_mtow = uconvert(u"m/s",InputValidate.get_value(df_mission,"Vmd_MTOW",col))
                    LD_max = InputValidate.get_value(df_mission,"LDmax",col)
                    V_md = Vmd_mtow.*sqrt.(α)

                    SFC_loiter = InputValidate.get_value(df_aircraft,"SFC_loiter",aircraft_idx)
                    SFC_cruise = InputValidate.get_value(df_aircraft,"SFC_cruise",aircraft_idx)

                    V = uconvert.(u"m/s", V_iter)
    
                    # SFC references
                    if engine_type != "Jet"
                        SFC_loiter = upreferred(SFC_loiter .* V)
                        SFC_cruise = upreferred(SFC_cruise .* V)
                    end

                    V_bar = V ./ V_md

                    R = uconvert(u"m", InputValidate.get_value(df_mission,"Distance",col))

                    range_save[ops_idx] = uconvert(u"km", R)

                    new_fraction =  exp.(upreferred.(-(R .* SFC_cruise .* g .* (V_bar .+ V_bar.^-3)) ./ (LD_max .* V_md .* 2)))
                    duration_save[:,save_idx] = uconvert.(u"minute", R ./ V)

                    # Save the old duration
                    duration_save_old[save_idx] = InputValidate.get_value(df_mission,"Duration",col)

                    save_idx += 1
                else
                    new_fraction = InputValidate.get_value(df_mission,"Fuel Fraction",col)
                end

                α = α .* new_fraction
            end

            # Get updated payload weight
            new_payload_weight = InputValidate.get_value(df_payload,"Payload Weight",payload_idx)

            # Obtain the fuel fraction
            new_fuel_frac = (1.0 .- α).*trapped_fuel

            # Calculate new fuel weight
            new_fuel_weight = (new_fuel_frac ./ (1 .- new_fuel_frac)) .* (new_payload_weight + empty_weight)
            
            if isnan(max_fuel_weight)
                continue
            end

            # See what velocities have fuel weights and payload below design one
            pass_fuel_idx = findall(value -> value <= max_fuel_weight, new_fuel_weight)
            pass_payload_idx = findall(value -> value <= max_payload_weight, new_payload_weight)

            # Ensure new fuel frac is not above 1 (or else crazy designs are allowed)
            pass_α_idx = findall(value -> value <= 1.0, new_fuel_frac)
            pass_fuel_idx = intersect(pass_fuel_idx,pass_α_idx)
            pass_payload_idx = intersect(pass_payload_idx,pass_α_idx)
            # If both fuel and payload is good, then calculate the minimum cost velocity!
            if !isempty(pass_fuel_idx) && !isempty(pass_payload_idx)
                pass_V = V_iter[pass_fuel_idx]
                pass_duration = duration_save[pass_fuel_idx,:]
                pass_fuel_weight = new_fuel_weight[pass_fuel_idx]
                cost_increase = 0.0
            else
                # Interpolate from payload-range frontier
                payload_cap = intp[range_save[ops_idx]]

                # How much less payload it could carry?
                payload_deficit = new_payload_weight - payload_cap

                # Get the velocity with the lowest deficit (i.e., not as much loss, it should be the minimum fuel speed!)
                pass_V = intp_V[range_save[ops_idx]]
                pass_duration = intp_duration[range_save[ops_idx]]
                pass_fuel_weight = max_fuel_weight

                # ASSUMPTION: Cost per kg per range per passenger not carried
                cost_per_pax_per_range = 0.15
                # This would also include range outside of the cruise range, which is fine, negligible
                total_range = sum(skipmissing(df_mission[findfirst(==("Distance"),df_mission[:,1]),ncol(df_mission)-N_stages+1:ncol(df_mission)]))
                pax_weight = InputValidate.get_value(df_payload,"Human Weight",payload_idx)

                cost_increase = cost_per_pax_per_range * (uconvert(u"kg",payload_deficit) / uconvert(u"kg",pax_weight)) * ustrip(uconvert(u"km",total_range))

                # Update weights
                df_payload = InputValidate.df_update_or_append(df=df_payload,label="Payload Weight",value=payload_cap,N_config=N_payload,col=payload_idx)
            end

            total_cost_opt = +Inf
            df_cost_breakdown_save = DataFrame()
            V_opt = NaN
            fuel_weight_opt = NaN

            # For each viable velocity
            for idx in eachindex(pass_V)
                save_idx = 1
                # Update mission profile to include new durations
                for col in key_cruise_monitor
                    duration = pass_duration[idx,save_idx]
                    df_mission = InputValidate.df_update_or_append(df=df_mission,label="Duration",value=duration,N_config=N_stages,col=col)

                    save_idx += 1
                end

                # Update with new fuel weight
                df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Fuel Weight",value=pass_fuel_weight[idx],N_config=N_payload,col=payload_idx)

                # Next, calculate the cost of operating and then return the averaged cost!
                (total_cost, df_cost_breakdown) = CostModel.cost_calc(df_cost_model = df_cost_model, df_aircraft = df_aircraft, aircraft_idx = aircraft_idx, df_mission = df_mission, N_stages = N_stages, df_payload = df_payload, payload_idx = payload_idx)

                # Add in the cost increase caused by not carrying some payload
                # Assumes ideal case where all passengers / payload could have been carried
                total_cost += cost_increase
                push!(df_cost_breakdown,["Cost of Missing Revenue","Lost Revenue",cost_increase])

                # Also add info about the payload and range
                push!(df_cost_breakdown,["Range","Range",cost_increase])
                push!(df_cost_breakdown,["Cost of Missing Revenue","Lost Revenue",cost_increase])

                # If minimum cost, then take the new value
                if total_cost < total_cost_opt
                    V_opt = pass_V[idx]
                    total_cost_opt = total_cost
                    df_cost_breakdown_save = df_cost_breakdown
                    fuel_weight_opt = pass_fuel_weight[idx]
                end
            end

            total_cost_save[ops_idx] = total_cost_opt
            V_optimal_save[ops_idx] = V_opt
            fuel_weight_save[ops_idx] = fuel_weight_opt

            if nrow(df_breakdown_save) == 0
                df_breakdown_save = df_cost_breakdown_save
                rename!(df_breakdown_save, "Cost (USD)" => "Scenario $ops_idx")
            else
                df_breakdown_save[:,"Scenario $ops_idx"] = df_cost_breakdown_save[:,"Cost (USD)"]
            end

            # Return to the old values
            df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Fuel Weight",value=max_fuel_weight,N_config=N_aircraft,col=aircraft_idx)
            df_payload = InputValidate.df_update_or_append(df=df_payload,label="Payload Weight",value=max_payload_weight,N_config=N_payload,col=payload_idx)

            save_idx = 1
            # Update mission profile to include old durations
            for col in key_cruise_monitor
                duration = duration_save_old[save_idx]
                df_mission = InputValidate.df_update_or_append(df=df_mission,label="Duration",value=duration,N_config=N_stages,col=col)

                save_idx += 1
            end
        end

        df_results = deepcopy(df_ops)
        df_results[:,"Payload"] = payload_save
        df_results[:,"Range"] = range_save
        df_results[:,"Total Cost"] = total_cost_save
        df_results[:,"Optimal Velocity"] = V_optimal_save
        df_results[:, "Fuel Weight"] = fuel_weight_save
        df_payloadrange = DataFrame(payload = payload_save, range = range_save)
        average_cost = mean(total_cost_save)

        return (average_cost, df_results, df_breakdown_save, save_payload_range)
    else
        # return nothing if no ops is given
        return (NaN, DataFrame(), DataFrame(), [])
    end
end

function plot_payload_range(df_save::DataFrame, row_idx::Int; metric = "Total")
    df_ops_results = df_save[row_idx, "Off-Mission Details"] 
    df_breakdown_save = df_save[row_idx, "Flight Cost Details"] 
    save_payload_range = df_save[row_idx, "Payload Range Diagram"]
    des_payload = df_save[row_idx, "Payload Weight"]
    des_range = df_save[row_idx, "Cruise Range"]

    p = Plots.plot(dpi = 500, legend = :bottomleft)

    # Get range, payload and costs
    range = try df_ops_results[:, "Range"]
    catch e
        @warn ("The iteration did not run successfully, cannot plot this!")
        return p
    end
    payload = df_ops_results[:, "Payload"]
    total_cost = df_ops_results[:, "Total Cost"]
    avg_total = round(mean(total_cost), sigdigits = 5)


    Plots.plot!(p, save_payload_range[1][1], save_payload_range[1][2], label = "Aircraft Performance Limit")
    Plots.scatter!(p, [des_range], [des_payload], label = "Design Mission, $des_range and $des_payload", markersize = 10, marker = :diamond, color = "black")
    Plots.scatter!(p, range, payload, zcolor = total_cost, c = cgrad(:RdYlGn, rev = true), clim = (10000, 80000), colorbar=true, colorbar_title="\n\nTotal Cost (USD)", label = "Off-Design Mission\nAverage = $avg_total USD",right_margin = 15Plots.mm)

    xlabel!("Range km")
    ylabel!("Payload kg")
    title!("Payload-Range Diagram")

    return p
end

end