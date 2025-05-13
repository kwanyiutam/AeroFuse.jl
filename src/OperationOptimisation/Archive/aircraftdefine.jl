# Functions that defines the aircraft properties
# e.g. Engine and Aircraft Configuration

##### Engine Definition
abstract type AbstractEngine end

mutable struct Engine{T} <: AbstractEngine
    type                :: String # type = Engine type, Jet or Propeller
    λ_bpr               :: T # λ_bpr = Engine bypass ratio (only for jet)
    prop_type           :: String # Propeller Type
    η_prop              :: T # η_prop = Propeller Efficiency (only for propeller)
    k_c                 :: T # k_c = coefficient based on fixed pitch or constant speed (only for propeller)
    prop_disc_load      :: T # prop_disc_load = Propeller disc loading (only for propeller), MUST be hp/ft^2
    SFC_cruise          :: T # SFC during cruise
    SFC_loiter          :: T # SFC during loiter
    T                   :: T # Thrust
end

function engine_set(type :: String, λ_bpr, η_prop, prop_type, prop_disc_load, SFC_cruise, SFC_loiter)
    k_c = NaN
    if type == "Jet"
        @assert λ_bpr > 0 "Bypass ratio must be positive"

        # TODO: Estimate SFC based on bypass ratio AND takeoff thrust estimation
        SFC_cruise = 14.1
        SFC_loiter = 11.3
    elseif type == "Propeller" || type == "Turboprop"
        @assert (η_prop >= 0 && η_prop <= 1) "Propeller efficiency must be between 0 and 1"

        if prop_type == "Constant Speed"
            k_c = 1
        elseif prop_type == "Fixed Pitch"
            k_c = 0.8
        else
            throw("Type of Propeller must be either constant speed or fixed pitch!")
        end

        SFC_cruise =  0.068
        SFC_loiter = 0.085
    else
        throw("Type of Engine must either be Jet or Propeller")
    end

    if type == "Turboprop"
        SFC_cruise =  0.085
        SFC_loiter = 0.101
    end

    T = promote_type(eltype(λ_bpr), eltype(η_prop), eltype(prop_disc_load), eltype(k_c), eltype(SFC_cruise), eltype(SFC_loiter))

    return(T, λ_bpr, η_prop, prop_disc_load, k_c, SFC_cruise, SFC_loiter)
end

function engine_param()
    return(["type", "λ_bpr", "η_prop", "prop_type", "prop_disc_load", "SFC_cruise", "SFC_loiter"])
end


function Engine(; type :: String = "Jet", λ_bpr = NaN, η_prop = NaN, prop_type :: String = NaN, prop_disc_load = NaN, SFC_cruise = NaN, SFC_loiter = NaN)

    (T, λ_bpr, η_prop, prop_disc_load, k_c, SFC_cruise, SFC_loiter) = engine_set(type, λ_bpr, η_prop, prop_type, prop_disc_load, SFC_cruise, SFC_loiter)

    return Engine{T}(type, λ_bpr, prop_type, η_prop, k_c, prop_disc_load, SFC_cruise, SFC_loiter, 0.0)
end

function Engine_Check(; engine :: Engine)

    (T, λ_bpr, η_prop, prop_disc_load, k_c, SFC_cruise, SFC_loiter) = engine_set(engine.type, engine.λ_bpr, engine.η_prop, engine.prop_type, engine.prop_disc_load, engine.SFC_cruise, engine.SFC_loiter)

    return Engine{T}(engine.type, λ_bpr, engine.prop_type, η_prop, k_c, prop_disc_load, SFC_cruise, SFC_loiter, 0.0)
end

##### Aircraft Definition
abstract type AbstractAircraft end

mutable struct Aircraft{T} <: AbstractAircraft
    engine      :: Engine{T}
    N_e         :: Int # N_e = Number of engines
    CL_max      :: T # CL_max = Maximum Coefficient of Lift
    ΔCL_max_to  :: T # ΔCL_max_to = Change in CL_max due to take-off flaps
    ΔCL_max_ld  :: T # ΔCL_max_ld = Change in CL_max due to landing flaps
    CD0         :: T # CD0 = Zero-lift drag
    ΔCD0_ldg    :: T # ΔCD0_ldg = Change in CD0 due to landing gear
    ΔCD0_tof    :: T # ΔCD0_tof = Change in CD0 due to take-off flaps
    ΔCD0_ldf    :: T # ΔCD0_ldf = Change in CD0 due to landing flaps
    AR          :: T # AR = Aspect ratio
    e           :: T # e = Oswald Efficiency
    Δe_ldg      :: T # Δe_ldg = Change in e due to landing gear
    Δe_tof      :: T # Δe_tof = Change in e due to take-off flaps
    Δe_ldf      :: T # Δe_ldf = Change in e due to landing flaps
    λ           :: T # λ = Taper ratio
    LD_max      :: T # LD_max = Maximum LD
    MTOW        :: T # Maximum takeoff weight
    Sref        :: T # Wing Area
    bref        :: T # Wing Span
    C_root      :: T # Root chord
    C_tip       :: T # Tip chord
end

function aircraft_set(;engine :: Engine, N_e :: Int, CL_max, ΔCL_max_to, ΔCL_max_ld, CD0, ΔCD0_ldg, ΔCD0_tof, ΔCD0_ldf, AR, e, Δe_ldg, Δe_tof, Δe_ldf, λ)
    T = promote_type(eltype(engine.λ_bpr), eltype(engine.η_prop), eltype(engine.prop_disc_load), eltype(engine.k_c),  eltype(N_e), eltype(CL_max), eltype(CD0), eltype(ΔCD0_ldg), eltype(ΔCD0_tof), eltype(ΔCD0_ldf), eltype(AR), eltype(e), eltype(Δe_ldg), eltype(Δe_tof), eltype(Δe_ldf), eltype(λ))
    @assert N_e > 0 "Number of engines must be bigger than 1."
    @assert CL_max > 0 && CD0 > 0 && AR > 0 "CL_max, CD0 and AR inputs must be positive"
    @assert e >= 0 && e <= 1 && λ >= 0 && λ <= 1 "e and λ must be between 0 and 1"
    
    LD_max = 0.5 * sqrt((pi * AR * e)/CD0)

    return(T, LD_max)
end

function aircraft_param()
    return(["N_e", "CL_max", "ΔCL_max_to", "ΔCL_max_ld", "CD0", "ΔCD0_ldg", "ΔCD0_tof", "ΔCD0_ldf", "AR", "e", "Δe_ldg", "Δe_tof", "Δe_ldf", "λ"])
end

function Aircraft(; engine :: Engine, N_e :: Int = 2, CL_max, ΔCL_max_to = 0.2, ΔCL_max_ld = 0.5, CD0, ΔCD0_ldg = 0.02, ΔCD0_tof = 0.02, ΔCD0_ldf = 0.07, AR, e, Δe_ldg = -0.05, Δe_tof = -0.05, Δe_ldf = -0.1, λ = 0.3)

    (T, LD_max) = aircraft_set(engine = engine, N_e = N_e,
    CL_max = CL_max, ΔCL_max_to = ΔCL_max_to, ΔCL_max_ld = ΔCL_max_ld, CD0 = CD0,
    ΔCD0_ldg = ΔCD0_ldg, ΔCD0_tof = ΔCD0_tof, ΔCD0_ldf = ΔCD0_ldf,
    AR = AR, e = e, Δe_ldg = Δe_ldg, Δe_tof = Δe_tof, Δe_ldf = Δe_ldf, λ = λ)
    
    return Aircraft{T}(engine, N_e, CL_max, ΔCL_max_to, ΔCL_max_ld, CD0, ΔCD0_ldg, ΔCD0_tof, ΔCD0_ldf, AR, e, Δe_ldg, Δe_tof, Δe_ldf, λ, LD_max, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function Aircraft_Check(; aircraft :: Aircraft)

    aircraft.engine = Engine_Check(engine = aircraft.engine)
    (T, LD_max) = aircraft_set(engine = aircraft.engine, N_e = aircraft.N_e,
    CL_max = aircraft.CL_max, ΔCL_max_to = aircraft.ΔCL_max_to, ΔCL_max_ld = aircraft.ΔCL_max_ld,
    CD0 = aircraft.CD0, ΔCD0_ldg = aircraft.ΔCD0_ldg, ΔCD0_tof = aircraft.ΔCD0_tof, ΔCD0_ldf = aircraft.ΔCD0_ldf,
    AR = aircraft.AR, e = aircraft.e, Δe_ldg = aircraft.Δe_ldg, Δe_tof = aircraft.Δe_tof, Δe_ldf = aircraft.Δe_ldf, λ = aircraft.λ)
    
    return Aircraft{T}(aircraft.engine, aircraft.N_e, aircraft.CL_max, aircraft.ΔCL_max_to, aircraft.ΔCL_max_ld, aircraft.CD0, aircraft.ΔCD0_ldg, aircraft.ΔCD0_tof, aircraft.ΔCD0_ldf, aircraft.AR, aircraft.e, aircraft.Δe_ldg, aircraft.Δe_tof, aircraft.Δe_ldf, aircraft.λ, LD_max, 0.0, 0.0, 0.0)
end

function change_aircraft(; aircraft :: Aircraft, config_input, idx_config, n_change)
    new_aircraft = deepcopy(aircraft)

    for idx_change in 1:n_change
        param_select = config_input[(idx_config-1)*2*n_change+(idx_change-1)*2+1]
        value_select = config_input[(idx_config-1)*2*n_change+idx_change*2]
        
        if param_select in engine_param()
            try
                setproperty!(new_aircraft.engine, Symbol(param_select), value_select)
            catch MethodError
                setproperty!(new_aircraft.engine, Symbol(param_select), parse(Float64, value_select))
            end
        else
            try
                setproperty!(new_aircraft, Symbol(param_select), value_select)
            catch MethodError
                setproperty!(new_aircraft, Symbol(param_select), parse(Float64, value_select))
            end
        end
    end

    new_aircraft = Aircraft_Check(aircraft = new_aircraft)

    return(new_aircraft)
end

function update_design(; aircraft :: Aircraft, WS, TW, MTOW)
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
    MTOW = MTOW * g
    aircraft.engine.T = ustrip(TW * MTOW)
    aircraft.MTOW = ustrip(MTOW)
    aircraft.Sref = ustrip(MTOW / WS)
    aircraft.bref = ustrip(sqrt(aircraft.AR * aircraft.Sref))
    aircraft.C_root = 2*aircraft.Sref/(aircraft.bref*(1+aircraft.λ)) # Using taper ratio to define 
    aircraft.C_tip = aircraft.C_root * aircraft.λ

    return aircraft
end