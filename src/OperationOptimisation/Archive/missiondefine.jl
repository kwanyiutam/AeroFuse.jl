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

function engine_init()
    engine_data = ParseData.get_data(data =  "JetEngine")

    SFC_data = engine_data[!, r"Bypass ratio|Cruise Sfc lb/hr/lb|Takeoff Thrust|Year"]
    SFC_data = dropmissing(SFC_data)

    SFC_data = filter(row -> row.Year > 1985, SFC_data)

    λ_bpr = SFC_data[:, 2]
    SFC = SFC_data[:, 3]
    T = SFC_data[:, 1]
    year = SFC_data[:, 4]
    SFC = SFC .* u"lb / hr / lbf"
    SFC = uconvert.(u"mg / N / s", SFC)

    SFC_nounit = ustrip(SFC)

    @. model(x, p) = p[1] + p[2] * exp(x[:,1] * p[3]) + p[4] * exp(x[:,2] * p[5])
    fit = curve_fit(model, hcat(T, λ_bpr), SFC_nounit, [120, 10, -0.001, -100, -0.1])

    return (model, fit, λ_bpr, SFC, T, year)
end

function range_weight_fractions(;engine :: Engine, range, range_unit :: String = "km", SFC, V, V_unit :: String = "knots", LD)
    ## TODO: Consider how LD depends on V

    @assert LD > 0 && V > 0 "L/D and V inputs must be positive"
    @assert SFC >= 0 && range >= 0 "SFC and range inputs must be non-negative"

    # Convert the input units into AeroUnits
    range_unit = AeroUnits.convert_to_unit(range_unit)
    V_unit = AeroUnits.convert_to_unit(V_unit)

    # Convert units
    range = uconvert(u"m", range * range_unit)
    V = uconvert(u"m/s", V * V_unit)
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant

    if engine.type == "Jet"
        SFC = SFC / 1000 * u"g / N / s" # Convert from mg to g
    else
        SFC = (SFC * u"mg / W / s") * V / engine.η_prop
    end

    W = exp((-range * SFC * g) / (V * LD))

    return W
end

function endurance_weight_fractions(;engine :: Engine, E, E_unit :: String = "minute", SFC, V, V_unit :: String = "knots", LD)
    ## TODO: Consider how LD depends on V

    @assert LD > 0 && V > 0 "L/D and V inputs must be positive"
    @assert SFC >= 0 && E >= 0 "SFC and endurance inputs must be non-negative"

    # Convert the input units into AeroUnits
    E_unit = AeroUnits.convert_to_unit(E_unit)
    V_unit = AeroUnits.convert_to_unit(V_unit)

    # Convert units
    E = uconvert(u"s", E * E_unit)
    V = uconvert(u"m/s", V * V_unit)
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant

    if engine.type == "Jet"
        SFC = SFC / 1000 * u"g / N / s" # Convert from mg to g
    else
        SFC = (SFC * u"mg / W / s") * V / engine.η_prop
    end

    W = exp((-E * SFC * g) / LD)

    return W
end

function takeoff_weight_fractions()
    # TODO: Model to predict weight fraction for takeoff
    return 0.97
end

function climb_weight_fractions()
    # TODO: Model to predict weight fraction for climb
    return 0.985
end

function descend_weight_fractions()
    # Optional for descend, quite insignificant
    return 0.995
end

##### Mission Profile Definition
abstract type AbstractMissionProfile end

mutable struct MissionProfile{T} <: AbstractMissionProfile
    aircraft    :: Aircraft{T}
    name        :: Array{String}
    stage       :: Array{String}
    w_frac      :: Array{T} # Weight fractions
    α           :: Array{T} # Weight fractions relative to MTOW
    h           :: Array{T} # Height
    V∞          :: Array{T} # Velocity
    distance    :: Array{T} # Distance travelled
    duration    :: Array{T} # Duration (endurance)
    G           :: Array{T} # Gradient (%)
    N_a         :: Array{Int64} # Engines available
end

function MissionProfile(; aircraft :: Aircraft)
    T = promote_type(eltype(aircraft.N_e), eltype(aircraft.CD0), eltype(aircraft.e), eltype(aircraft.engine.η_prop), eltype(aircraft.engine.λ_bpr))
    
    name = []
    stage = []
    w_frac = []
    α = []
    h = []
    V∞ = []
    distance = []
    duration = []
    G = []
    N_a = []

    return MissionProfile{T}(aircraft, name, stage, w_frac, α, h, V∞, distance, duration, G, N_a)
end

function add_mission(; mission :: MissionProfile, name :: String, add_stage :: String, order = 0, V = 0.0, V_unit :: String = "knots", h = 0.0, h_unit :: String = "ft", distance = 0.0, distance_unit :: String = "km", duration = 0.0, duration_unit :: String = "minute", G = 0.0, N_a :: Int64 = 2)
    current_stage = mission.stage

    @assert order >= 0 && order <= (length(current_stage) + 1) "Order must be 0 (append stage) or smaller than or equal to the number of existing stages + 1"
    @assert distance >= 0 && duration >= 0 && V >= 0 && h >= 0 "Distance, Duration, Heights and Velocities must be bigger than zero"
    @assert N_a >= 1 && N_a <= mission.aircraft.N_e "Number of available engines must be bigger than 1 and less than or equal to number of engines."

    # Convert units into Unitful units
    h_unit = AeroUnits.convert_to_unit(h_unit)
    V_unit = AeroUnits.convert_to_unit(V_unit)
    distance_unit = AeroUnits.convert_to_unit(distance_unit)
    duration_unit = AeroUnits.convert_to_unit(duration_unit)

    h = uconvert(u"m", h * h_unit)
    V = uconvert(u"m/s", V * V_unit)
    distance = uconvert(u"m", distance * distance_unit)
    duration = uconvert(u"s", duration * duration_unit)

    if order == 0
        order = length(current_stage) + 1
    end

    if add_stage == "Takeoff"
        check_to = findall(y -> y == order, findall(x -> x == "Descend", current_stage).+ 1)
        @assert order == 1 || length(check_to) > 0 "Order of takeoff stage must be either the first stage or after a descend phase"

        if duration == 0*u"s"
            duration = uconvert(u"s", 15 * u"minute") # Estimate duration for taxi + takeoff time
        end

        w_frac = takeoff_weight_fractions()
    elseif add_stage == "Climb"
        check_cl = findall(y -> y == order, findall(x -> x == "Takeoff" || x == "Descend", current_stage).+ 1)
        @assert length(check_cl) > 0 "Order of climb stage must be after a takeoff or descend phase"
        @assert G > 0 "Climb must have a positive gradient"

        if duration == 0*u"s"
            duration = uconvert(u"s", 15 * u"minute") # Estimate duration for climb time
        end
        
        w_frac = climb_weight_fractions()
    elseif add_stage == "Cruise"
        check_cr = findall(y -> y == order, (findall(x -> x == "Climb" || x == "Descend" || x == "Loiter", current_stage).+ 1))
        @assert length(check_cr) > 0 "Order of cruise stage must be after a climb, descend or loiter phase"
        @assert G == 0 "Cruise must have a zero gradient"
        @assert distance > 0*u"m" "Distance must be longer than zero"

        if mission.aircraft.engine.type == "Jet"
            LD = mission.aircraft.LD_max * 0.866
        else
            LD = mission.aircraft.LD_max
        end

        duration = uconvert(u"s", distance / V) #  Override duration time
        w_frac = range_weight_fractions(engine = mission.aircraft.engine, range = ustrip(distance), range_unit = "m", SFC = mission.aircraft.engine.SFC_cruise, V = ustrip(V), V_unit = "m/s", LD = LD)
    elseif add_stage == "Descend"
        check_ds = findall(y -> y == order, findall(x -> x == "Climb" || x == "Cruise" || x == "Loiter", current_stage).+ 1)
        @assert length(check_ds) > 0 "Order of descend stage must be after a climb, cruise or loiter phase"
        @assert G < 0 "Descend must have a negative gradient"

        if duration == 0*u"s"
            duration = uconvert(u"s", 20 * u"minute") # Estimate duration for descend time
        end
        
        w_frac = descend_weight_fractions()
    elseif add_stage == "Loiter"
        check_lo = findall(y -> y == order, (findall(x -> x == "Climb" || x == "Descend" || x == "Cruise" , current_stage).+ 1))
        @assert length(check_lo) > 0 "Order of loiter stage must be after a climb, descend or cruise phase"
        @assert G == 0 "Loiter must have a zero gradient"
        @assert duration > 0*u"s" "Duration must be longer than zero"

        if mission.aircraft.engine.type == "Jet"
            LD = mission.aircraft.LD_max
        else
            LD = mission.aircraft.LD_max * 0.866
        end

        w_frac = endurance_weight_fractions(engine = mission.aircraft.engine, E = ustrip(duration), E_unit = "s", SFC = mission.aircraft.engine.SFC_loiter, V = ustrip(V), V_unit = "m/s", LD = LD)
    elseif add_stage == "Landing"
        check_ld = findall(y -> y == order, (findall(x -> x == "Descend", current_stage).+ 1))
        @assert length(check_ld) > 0 "Order of landing stage must be after a descend phase"

        if duration == 0*u"s"
            duration = duration = uconvert(u"s", 15 * u"minute") # Estimate duration for landing time
        end

        w_frac = 1
    else
        # Return error if the additional stage is not valid
        throw("Additional stage is not valid")
    end

    # Strip units
    h = ustrip(h)
    V = ustrip(V)
    distance = ustrip(distance)
    duration = ustrip(duration)


    # Mutate the mission profile struct
    mission.name = insert!(mission.name, order, name)
    mission.stage = insert!(mission.stage, order, add_stage)
    mission.h = insert!(mission.h, order, h)
    mission.V∞ = insert!(mission.V∞, order, V)
    mission.distance = insert!(mission.distance, order, distance)
    mission.duration = insert!(mission.duration, order, duration)
    mission.w_frac = insert!(mission.w_frac, order, w_frac)
    mission.G = insert!(mission.G, order, G)
    mission.N_a = insert!(mission.N_a, order, N_a)

    if order == 1
        mission.α = insert!(mission.α, order, w_frac)
    else
        mission.α = insert!(mission.α, order, mission.α[order - 1]*w_frac)
    end

    return mission
end

##### Payload Definition
abstract type AbstractPayload end

mutable struct Payload{T} <: AbstractPayload
    N_pax           :: Int64 # Passengers
    N_fcrew         :: Int64 # Flight Crew (Should be based on regulation later)
    N_ccrew         :: Int64 # Cabin Crew (Should be based on regulation later)
    W_drop          :: T     # Droppable Payload
    W_payload       :: T     # Total payload weight
end

function Payload(; N_pax :: Int64 = 0, N_fcrew :: Int64 = 2, N_ccrew :: Int64 = 4, W_human = 80.0, W_human_unit :: String = "kg", W_bag = 15.0, W_bag_unit :: String = "kg", W_drop = 0.0, W_drop_unit :: String = "kg", W_other = 0.0, W_other_unit :: String = "kg")
    @assert N_pax >= 0 && N_fcrew >= 0 && N_ccrew >= 0 "Number of people must be non-negative!"
    @assert W_human >= 0 && W_bag >= 0 && W_drop >= 0 && W_other >= 0 "Weight must be non-negative!"
    
    T = promote_type(eltype(W_other), eltype(W_bag), eltype(W_drop))

    # Convert units into Unitful units
    W_human_unit = AeroUnits.convert_to_unit(W_human_unit)
    W_bag_unit = AeroUnits.convert_to_unit(W_bag_unit)
    W_drop_unit = AeroUnits.convert_to_unit(W_drop_unit)
    W_other_unit = AeroUnits.convert_to_unit(W_other_unit)
    
    W_drop = ustrip(uconvert(u"kg", W_drop * W_drop_unit))

    N_human = N_pax + N_fcrew + N_ccrew
    W_payload = N_human * (ustrip(uconvert(u"kg", W_human * W_human_unit)) + ustrip(uconvert(u"kg", W_bag * W_bag_unit))) + ustrip(uconvert(u"kg", W_other * W_other_unit)) + W_drop

    return Payload{T}(N_pax, N_fcrew, N_ccrew, W_drop, W_payload)
end
