# Initialise packages used
using CSV
using DataFrames
using LsqFit
using Unitful
using ISAData

using Plots

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the ParseData module
include("parsedata.jl")

# Include the AircraftDefine functions
include("aircraftdefine.jl")

function parse_DOC()
    DOC_data = ParseData.get_data(data = "DOC")
    inflation_data = ParseData.get_data(data =  "Inflation")

    master_data = leftjoin(DOC_data[1:nrow(DOC_data) .!= only(findall(==("Units"), DOC_data.Category)),:], inflation_data, on = :Year)

    predictors = names(DOC_data[!, r"^((?!(Data|Category|Year|Constant)).)*$"])
    units = Array(DOC_data[only(findall(==("Units"), DOC_data.Category)), r"^((?!(Data|Category|Year|Constant)).)*$"])
    predictors = Dict(predictors .=> units)
    parameters = (master_data[!, r"^(Data|Category)"])
    
    return (parameters, predictors, master_data, inflation_data)
end

function crew_cost(;payload :: Payload, W_0, W_0_unit :: String = "kg", time = 1, time_unit :: String = "hr", inflation)
    # Get inflation rate compared to latest figures
    year_data = 2010 # Data origin's year
    factor = inflation[nrow(inflation),"CPI"]./inflation[inflation.Year .== year_data,"CPI"] # Last row = latest value

    # Convert units into Unitful units
    W_0_unit = AeroUnits.convert_to_unit(W_0_unit)
    time_unit = AeroUnits.convert_to_unit(time_unit)

    # Convert to correct units
    W_0 = ustrip(uconvert.(u"kg", W_0 * W_0_unit))
    time = ustrip(uconvert.(u"hr", time * time_unit))

    # Calculate the costs, including inflation to current year
    cost_cabin = 81 .* time .* payload.N_ccrew .* factor
    cost_flight = (0.326*W_0/1000 .+ 653) .* time .* factor

    return (cost_cabin, cost_flight)
end

function depreciation_cost(;N_e, W_0, W_0_unit :: String = "kg", W_e, W_e_unit :: String = "kg", T, T_unit :: String = "N", time = 1, time_unit :: String = "hr", inflation)
    # Get inflation rate compared to latest figures
    year_data = 2010 # Data origin's year
    factor = inflation[nrow(inflation),"CPI"]./inflation[inflation.Year .== year_data,"CPI"] # Last row = latest value

    # Convert units into Unitful units
    W_0_unit = AeroUnits.convert_to_unit(W_0_unit)
    W_e_unit = AeroUnits.convert_to_unit(W_e_unit)
    T_unit = AeroUnits.convert_to_unit(T_unit)
    time_unit = AeroUnits.convert_to_unit(time_unit)

    # Convert to correct units
    W_0 = ustrip.(uconvert.(u"kg", W_0 * W_0_unit))
    W_e = ustrip.(uconvert.(u"kg", W_e * W_e_unit))
    T = ustrip.(uconvert.(u"lbf", T * T_unit))
    time = ustrip.(uconvert.(u"hr", time * time_unit))

    utilisation = 6100 .- (3100 .*time.^(-0.3342))

    ac_cost = zeros(length(W_e), 1)

    for idx in 1:length(W_e)
        if W_e[idx,1] >= 10000
            ac_cost[idx] = 10^6 * (1.18 * W_e[idx]^0.48 - 116)
        else
            ac_cost[idx] = -0.002695 * W_e[idx]^2 + 1967 * W_e[idx] - 2158000
        end
    end

    ac_cost = ac_cost .* factor
    eng_cost = 1.76 .* 82.5 .* T .* N_e .* factor # TODO: What about turboprop?
    af_cost = ac_cost .- eng_cost

    deprec_cost = (0.9 .* time).*(ac_cost .+ 0.1 .* af_cost .+ 0.3 .* eng_cost) ./ (14 .* utilisation)

    return(deprec_cost, ac_cost, eng_cost, af_cost, utilisation)
end

function fuel_calculate(;aircraft:: Aircraft, V_md, V_md_unit :: String = "m/s", V_cruise, V_cruise_unit :: String = "m/s", cruise_range, cruise_range_unit :: String = "km")
    # Convert units into Unitful units
    V_md_unit = AeroUnits.convert_to_unit(V_md_unit)
    V_cruise_unit = AeroUnits.convert_to_unit(V_cruise_unit)
    cruise_range_unit = AeroUnits.convert_to_unit(cruise_range_unit)

    V_md = uconvert(u"m/s", V_md * V_md_unit)
    V_cruise = uconvert.(u"m/s", V_cruise * V_cruise_unit)
    cruise_range = uconvert(u"m", cruise_range * cruise_range_unit)

    V_bar = V_cruise ./ V_md

    SFC = aircraft.engine.SFC_cruise

    if aircraft.engine.type == "Jet"
        SFC = SFC / 1000 * u"g / N / s" # Convert from mg to g
    else
        SFC = (SFC / 1000 * u"g / W / s") .* V_cruise ./ aircraft.engine.η_prop
    end

    # Calculate the fuel costs
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
    fuel_weight_ratio = exp.(upreferred.((-cruise_range * SFC * g .* (V_bar .+ (V_bar.^(-3)))) ./ (V_md * aircraft.LD_max * 2)))
    
    return fuel_weight_ratio
end

function velocity_choose(;design_input :: String, aircraft:: Aircraft, V_md, V_md_unit :: String = "m/s", V_user, V_user_unit :: String = "m/s")
    # Convert units into Unitful units
    V_md_unit = AeroUnits.convert_to_unit(V_md_unit)
    V_user_unit = AeroUnits.convert_to_unit(V_user_unit)

    V_md = uconvert(u"m/s", V_md * V_md_unit)
    V_user = uconvert(u"m/s", V_user * V_user_unit)

    if design_input == "Maximum Range" || design_input == "Minimum Cost"
        if aircraft.engine.type == "Jet"
            V_cruise = uconvert.(u"knots", V_md * 3^(1/4))
        else
            V_cruise = uconvert.(u"knots", V_md)
        end
    elseif design_input == "User-Define"
        V_cruise = uconvert.(u"knots", V_user)
    else
        error("Invalid Design Input Value")
    end

    V_convert = uconvert.(u"m/s", V_cruise)

    return V_convert
end

function cost_calculate(;aircraft :: Aircraft, payload :: Payload, cruise_range, cruise_range_unit :: String = "km", V_cruise, V_cruise_unit :: String = "m/s", V_md, V_md_unit :: String = "m/s", W_f_other, W_0, W_0_unit :: String = "kg", W_e, W_e_unit :: String = "kg", time_others, inflation)
    W_0_unit = AeroUnits.convert_to_unit(W_0_unit)
    W_e_unit = AeroUnits.convert_to_unit(W_e_unit)
    cruise_range_unit = AeroUnits.convert_to_unit(cruise_range_unit)
    V_cruise_unit = AeroUnits.convert_to_unit(V_cruise_unit)
    V_md_unit = AeroUnits.convert_to_unit(V_md_unit)

    W_0 = uconvert.(u"kg", W_0 * W_0_unit)
    W_e = uconvert.(u"kg", W_e * W_e_unit)
    cruise_range = uconvert(u"m", cruise_range * cruise_range_unit)
    V_cruise = uconvert.(u"m/s", V_cruise * V_cruise_unit)
    V_md = uconvert(u"m/s", V_md * V_md_unit)

    # Assuming takeoff, landing, taxi is constant regardless of aircraft design, only care about cruise
    t_f = cruise_range ./ V_cruise
    
    fuel_weight_ratio = fuel_calculate(
        aircraft = aircraft,
        V_md = ustrip(V_md),
        V_cruise = ustrip(V_cruise),
        cruise_range = ustrip(cruise_range),
        cruise_range_unit = "m")
    
    W_f_total = ((1 .- fuel_weight_ratio) .+ (1 .- W_f_other))
    fuel_weight = W_f_total .* W_0 # multiply new fuel weight ratio
    fuel_cost = fuel_weight * (691.95/1000 * u"kg^-1") # Fuel cost, not always valid 691.95

    # Calculate the crew costs
    (cost_ccrew, cost_fcrew) = crew_cost(payload = payload, W_0 = ustrip(W_0), time = ustrip(t_f) .+ time_others, time_unit = "s", inflation = inflation)
    
    # Calculate deprecaition cost
    (cost_depr, _, _, _, _) = depreciation_cost(N_e = aircraft.N_e,W_0 = ustrip(W_0), W_e = ustrip(W_e), T = 0, time = ustrip(t_f) .+ time_others, time_unit = "s", inflation = inflation)
    
    # Calculate the maintenance costs
    # TODO

    total_cost = cost_ccrew .+ cost_fcrew .+ fuel_cost .+ cost_depr

    return(t_f, total_cost, fuel_cost, cost_ccrew, cost_fcrew, cost_depr, fuel_weight_ratio)
end

function weight_converge(;W_0, W_payload, new_fuel_final, A, C)
    W_e = NaN
    W_0_prev = 0
    max_iter = 50
    iter = 1

    while abs(W_0-W_0_prev) > 0.01
        W_e = A*(W_0^C) # Empty weight Estimate

        W_0_prev = W_0
        W_0 = (W_payload) / (1 - new_fuel_final - W_e)

        iter = iter + 1

        if W_0 < 0
            W_0 = NaN
            W_e = NaN
            break
        elseif iter >= max_iter
            break
        end
    end

    return(W_0, W_e)
end

function weight_estimate(;
    stage_idx :: Int64,
    design_input :: String,
    mission :: MissionProfile,
    payload :: Payload,
    WS_min = 1,
    WS_max = 10000,
    WS_N = 5000,
    T_min = 1,
    T_max = 1000000,
    T_max_N = 1000,
    α_correction,
    time_correction,
    A,
    C,
    W_0_init,
    DOC_data,
    inflation
    )
    
    aircraft = mission.aircraft
    h = mission.h[stage_idx]
    α = mission.α[stage_idx - 1] .* α_correction # Get the previous weight fraction for maximum alpha
    V_user = mission.V∞[stage_idx]
    old_W_frac = mission.w_frac[stage_idx]
    old_α_final = mission.α[length(mission.α)] .* α_correction
    α_others = old_α_final ./ old_W_frac # Other fuel weights
    time_others = sum(mission.duration[1:end .!= stage_idx]) .+ time_correction
    cruise_range = mission.distance[stage_idx]
    fuel_trap = 1.06 # Trapped fuel multiplier

    (ρ_0,_,_,_) = ISAdata(0*u"m")
    (ρ,_,_,_) = ISAdata(h)
    σ = ρ / ρ_0

    WS = LinRange(WS_min*u"kg/m/s^2", WS_max*u"kg/m/s^2", WS_N)

    V_imd = sqrt.(2*α.*WS./ρ_0).*((1/(pi*aircraft.AR*aircraft.e*aircraft.CD0)).^0.25)
    V_md = V_imd ./ sqrt(σ)

    W_0_save = zeros(length(WS),1)
    W_e_save = zeros(length(WS),1)
    W_frac_save = zeros(length(WS),1)
    W_f_save = zeros(length(WS),1)
    V_save = zeros(length(WS),1)
    total_cost_save = zeros(length(WS),T_max_N)
    cost_ccrew_save = zeros(length(WS),T_max_N)
    cost_fcrew_save = zeros(length(WS),T_max_N)
    cost_depr_save = zeros(length(WS),T_max_N)
    fuel_cost_save = zeros(length(WS),T_max_N)
    time_elapsed_save = zeros(length(WS),T_max_N)
    W_payload = payload.W_payload

    for idx in 1:length(WS)
        V_md_idx = ustrip(V_md[idx])

        if length(α_others) > 1
            α_other = α_others[idx,1]
            time_other = time_others[idx,1]
            old_α_fin = old_α_final[idx,1]
        else
            α_other = α_others[1]
            time_other = time_others[1]
            old_α_fin = old_α_final[1]
        end

        if isnan(α_other)
            W_0_save[idx] = NaN
            W_e_save[idx] = NaN
            W_frac_save[idx] = NaN
            W_f_save[idx] = NaN
            V_save[idx] = NaN
            total_cost_save[idx,:] .= NaN
            cost_ccrew_save[idx] = NaN
            cost_fcrew_save[idx] = NaN
            cost_depr_save[idx,:] .= NaN
            fuel_cost_save[idx] = NaN
            time_elapsed_save[idx] = NaN
            continue
        end

        if design_input == "Minimum Cost"
            global W_0 = NaN
            global W_e = NaN
            global W_f = NaN
            global W_frac = NaN
            global V_min_cost = NaN
            global total_cost_new = 9*10^20
            global cost_ccrew_new = NaN
            global cost_fcrew_new = NaN
            global cost_depr_new = NaN
            global fuel_cost_new = NaN
            global time_elapse_new = NaN

            V_iterate = LinRange(1*u"knots", 550*u"knots", 999)
            V_convert = ustrip.(uconvert.(u"m/s", V_iterate))

            fuel_w_frac = fuel_calculate(
                aircraft = aircraft,
                V_md = V_md_idx,
                V_cruise = V_convert,
                cruise_range = cruise_range,
                cruise_range_unit = "m"
            )

            W_0_save_2 = zeros(length(V_iterate),1)
            W_e_save_2 = zeros(length(V_iterate),1)
            new_fuel_final_save = zeros(length(V_iterate),1)

            for j in 1:lastindex(V_convert)
                new_fuel_final_save[j] = (1 - (old_α_fin / old_W_frac * fuel_w_frac[j])) * fuel_trap # Update the alpha

                (W_0_save_2[j], W_e_save_2[j]) = weight_converge(W_0 = W_0_init, W_payload = W_payload, new_fuel_final = new_fuel_final_save[j], A = A, C = C)
            end

            W_e_save_2 = W_e_save_2 .* W_0_save_2

            (time_elapse, total_cost, fuel_cost, cost_ccrew, cost_fcrew, cost_depr, _) = cost_calculate(;
                aircraft = aircraft,
                payload = payload,
                cruise_range = cruise_range,
                cruise_range_unit = "m",
                V_cruise = V_convert,
                V_cruise_unit = "m/s",
                V_md = V_md_idx,
                W_f_other = α_other,
                W_0 = W_0_save_2,
                W_e = W_e_save_2,
                time_others = time_other,
                inflation = inflation
            )

            total_cost[isnan.(total_cost)] .= +Inf
            idx_opt = argmin(total_cost)
            total_cost[total_cost .== +Inf] .= NaN

            # If all values were NaN, then just give up
            if sum(isnan.(total_cost)) == length(total_cost)
                global W_frac = NaN
                global W_0 = NaN
                global W_e = NaN
                global W_f = NaN
                global V_min_cost = NaN
                global fuel_cost_new = NaN
                global cost_ccrew_new = NaN
                global cost_fcrew_new = NaN
                global cost_depr_new = fill(NaN, T_max_N)
                global time_elapse_new = NaN
                global total_cost_new = fill(NaN, T_max_N)
            else
                global W_frac = fuel_w_frac[idx_opt, 1]
                global W_0 = W_0_save_2[idx_opt, 1]
                global W_e = W_e_save_2[idx_opt, 1]
                global W_f = new_fuel_final_save[idx_opt, 1]
                global V_min_cost = V_iterate[idx_opt, 1]
                global fuel_cost_new = fuel_cost[idx_opt, 1]
                global cost_ccrew_new = cost_ccrew[idx_opt, 1]
                global cost_fcrew_new = cost_fcrew[idx_opt, 1]
                global time_elapse_new = time_elapse[idx_opt, 1]

                T_range = LinRange(0, 0, T_max_N) # TODO:  LinRange(T_min, T_max, T_max_N) TEMP change to 0, so that propeller doesn't always lose
                (cost_depr_update, _, _, _, _) = depreciation_cost(N_e = aircraft.N_e, W_0 = ustrip(W_0), W_e = ustrip(W_e), T = T_range, time = ustrip(time_elapse_new) + time_other, time_unit = "s", inflation = inflation)

                global cost_depr_new = cost_depr_update
                global total_cost_new = total_cost[idx_opt, 1] .- cost_depr[idx_opt, 1] .+ cost_depr_update
            end
        else
            # Skip the stuff if not minimising cost
            V_cruise = velocity_choose(
                design_input = design_input,
                aircraft = aircraft,
                V_md = V_md_idx,
                V_user = V_user
            )

            fuel_w_frac = fuel_calculate(
                aircraft = aircraft,
                V_md = V_md_idx,
                V_cruise = ustrip(V_cruise),
                cruise_range = cruise_range,
                cruise_range_unit = "m"
            )

            new_fuel_final = (1 - (old_α_fin / old_W_frac * fuel_w_frac)) * fuel_trap # Update the alpha

            (W_0, W_e) = weight_converge(W_0 = W_0_init, W_payload = W_payload, new_fuel_final = new_fuel_final, A = A, C = C)
            
            W_e = W_e * W_0

            (time_elapse, total_cost, fuel_cost, cost_ccrew, cost_fcrew, cost_depr, _) = cost_calculate(;
                aircraft = aircraft,
                payload = payload,
                cruise_range = cruise_range,
                cruise_range_unit = "m",
                V_cruise = ustrip(V_cruise),
                V_cruise_unit = "m/s",
                V_md = V_md_idx,
                W_f_other = α_other,
                W_0 = W_0,
                W_e = W_e,
                time_others = time_other,
                inflation = inflation
            )

            global W_frac = fuel_w_frac
            global W_0 = W_0
            global W_e = W_e
            global W_f = new_fuel_final
            global V_min_cost = uconvert(u"knots", V_cruise)
            global fuel_cost_new = fuel_cost[1]
            global cost_ccrew_new = cost_ccrew[1]
            global cost_fcrew_new = cost_fcrew[1]
            global time_elapse_new = time_elapse[1]

            T_range = LinRange(0, 0, T_max_N) # TODO:  LinRange(T_min, T_max, T_max_N) TEMP change to 0, so that propeller doesn't always lose
            (cost_depr_update, _, _, _, _) = depreciation_cost(N_e = aircraft.N_e, W_0 = ustrip(W_0), W_e = ustrip(W_e), T = T_range, time = ustrip(time_elapse_new) .+ time_other, time_unit = "s", inflation = inflation)

            global cost_depr_new = cost_depr_update
            global total_cost_new = total_cost[1] .- cost_depr[1] .+ cost_depr_update
        end

        W_0_save[idx] = W_0
        W_e_save[idx] = W_e
        W_frac_save[idx] = W_frac
        W_f_save[idx] = W_f
        V_save[idx] = ustrip(V_min_cost)
        total_cost_save[idx,:] = total_cost_new
        cost_ccrew_save[idx,:] .= cost_ccrew_new
        cost_fcrew_save[idx,:] .= cost_fcrew_new
        cost_depr_save[idx,:] = cost_depr_new
        fuel_cost_save[idx,:] .= fuel_cost_new
        time_elapsed_save[idx,:] .= ustrip(time_elapse_new)
    end

    α_correction = α_correction .* (W_frac_save ./ old_W_frac) # Final alpha
    time_correction = time_correction .+ time_elapsed_save .- mission.duration[stage_idx]

    weight_result = (
        W_0 = W_0_save,
        W_e_frac = W_e_save ./ W_0_save,
        W_f_frac = W_f_save
    )

    cost_result = (
        total_cost = total_cost_save,
        ccrew_cost = cost_ccrew_save,
        fcrew_cost = cost_fcrew_save,
        depr_cost = cost_depr_save,
        fuel_cost = fuel_cost_save
    )

    return (weight_result, cost_result, V_save, W_frac_save, α_correction, time_correction)
end
