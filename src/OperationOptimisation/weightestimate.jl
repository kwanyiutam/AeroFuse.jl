module WeightEst

# Initialise packages used
using DataFrames
using Unitful
using AeroFuse

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the InputValidate module
include("inputvalidate.jl")

# Include AircraftAero module
include("aircraftaero.jl")

function wing_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    # Key parameters for weight estimation
    W_dg = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx) # Assume same as MTOW
    N_z = InputValidate.get_value(df_aircraft,"Load Factor Positive",aircraft_idx) * 1.5 # Assume same as limit load factor x 1.5
    S_W = InputValidate.get_value(df_aircraft,"Wing Area",aircraft_idx) # Wing Sref
    AR = InputValidate.get_value(df_aircraft,"Wing AR",aircraft_idx) # Wing AR
    λ = InputValidate.get_value(df_aircraft,"Wing Taper Ratio",aircraft_idx) # Wing AR
    b_W = InputValidate.get_value(df_aircraft,"Wing Span",aircraft_idx) # Wing AR

    # Calculate control surface area
    c_root = AircraftAero.root_chord(S_W,b_W,λ)
    y_half_1 = 0.1 # ASSUMPTION! Control surface starts at 10% half span
    y_half_2 = 0.9 # ASSUMPTION! Control surface ends at 90% half span
    c_1 = AircraftAero.chord_pos(y_half_1,λ,c_root)
    c_2 = AircraftAero.chord_pos(y_half_2,λ,c_root)
    control_chord_proportion = 0.3 # ASSUMPTION! control surface takes up 30% of wing chord
    S_csw = (c_1+c_2)*control_chord_proportion/2*(y_half_2-y_half_1)*(b_W/2)
    
    # Weight estimation
    Λ = InputValidate.get_value(df_aircraft,"Wing Quarterchord Sweep",aircraft_idx) # Wing Quarterchord Sweep
    foil = AircraftAero.parse_aerofoil(InputValidate.get_value(df_aircraft,"Wing Airfoil",aircraft_idx))
    t_c = AeroFuse.maximum_thickness_to_chord(foil)[2]

    # Conversion to imperial
    W_dg = uconvert(u"lb", W_dg)
    S_W = uconvert(u"ft^2", S_W)
    b_W = uconvert(u"ft", b_W)
    S_csw = uconvert(u"ft^2", S_csw)

    # Calculate the weight
    W_wing = 0.0051*((W_dg*N_z)^0.557*S_W^0.649*AR^0.5*(1+λ)^0.1*S_csw^0.1)
    W_wing /= (cosd(Λ) * t_c^0.4)

    return uconvert(u"kg", ustrip(W_wing) * u"lb")
end

function horizontal_tail_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    # Get required parameters from data
    K_uht = 1.0 # Assumption: not all moving tail
    W_dg = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx)
    N_z = InputValidate.get_value(df_aircraft,"Load Factor Positive",aircraft_idx) * 1.5
    S_ht = InputValidate.get_value(df_aircraft,"HT Area",aircraft_idx)  # Horizontal tail area
    AR_h = InputValidate.get_value(df_aircraft,"HT AR",aircraft_idx)    # HT aspect ratio
    Λ_ht = InputValidate.get_value(df_aircraft,"HT Quarterchord Sweep",aircraft_idx) # degrees
    L_ht = InputValidate.get_value(df_aircraft,"HT Tail Arm",aircraft_idx)  # Tail moment arm
    F_w = InputValidate.get_value(df_aircraft,"Diameter",aircraft_idx) / 10 # ASSUMPTION! the fuselage width is 1/10 of the diameter
    B_h = InputValidate.get_value(df_aircraft,"HT Span",aircraft_idx)
    K_y = 0.3 * L_ht  # ASSUMPTION! radius of gyration constant (typical values: 0.3)
    S_e = 0.2 * S_ht  # ASSUMPTION! elevator area ~20% of tail area

    # Convert to imperial
    W_dg = uconvert(u"lb", W_dg)
    S_ht = uconvert(u"ft^2", S_ht)
    L_ht = uconvert(u"ft", L_ht)
    F_w = uconvert(u"ft", F_w)
    B_h = uconvert(u"ft", B_h)
    K_y = uconvert(u"ft", K_y)
    S_e = uconvert(u"ft^2", S_e)

    # Compute the weight
    W_ht = 0.0379 * K_uht * W_dg^0.639 * N_z^0.1 * S_ht^0.75 * K_y^0.704 * AR_h^0.166 * (1 + S_e/S_ht)^0.1
    W_ht /= ((1 + F_w/B_h)^0.25 * L_ht * cosd(Λ_ht))

    return uconvert(u"kg", ustrip(W_ht) * u"lb")
end

function vertical_tail_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    # Input values
    W_dg = InputValidate.get_value(df_aircraft, "MTOW", aircraft_idx)
    N_z = InputValidate.get_value(df_aircraft, "Load Factor Positive", aircraft_idx) * 1.5
    S_vt = InputValidate.get_value(df_aircraft, "VT Area", aircraft_idx)
    AR_v = InputValidate.get_value(df_aircraft, "VT AR", aircraft_idx)
    Λ_vt = InputValidate.get_value(df_aircraft, "VT Quarterchord Sweep", aircraft_idx)
    L_vt = InputValidate.get_value(df_aircraft, "VT Tail Arm", aircraft_idx)
    H_t_H_v = InputValidate.get_value(df_aircraft, "HT Root Chord Location (Relative to Span of VT)", aircraft_idx) # Vertical location along the span of VT
    K_z = L_vt   # ASSUMPTION! Equal to the VT tail arm

    # Obtain thickness to chord ratio
    foil = AircraftAero.parse_aerofoil(InputValidate.get_value(df_aircraft,"VT Airfoil",aircraft_idx))
    t_c_v = AeroFuse.maximum_thickness_to_chord(foil)[2]

    # Unit conversions
    W_dg = uconvert(u"lb", W_dg)
    S_vt = uconvert(u"ft^2", S_vt)
    L_vt = uconvert(u"ft", L_vt)
    K_z = uconvert(u"ft", K_z)

    # Weight calculation
    W_vt = 0.0026 * (1 + H_t_H_v)^0.225 * W_dg^0.556 * N_z^0.536 * S_vt^0.5 * K_z^0.875 * AR_v^0.35
    W_vt /= (L_vt^0.5 * cosd(Λ_vt) * t_c_v^0.5)

    return uconvert(u"kg", ustrip(W_vt) * u"lb")
end

function fuselage_weight(df_aircraft, aircraft_idx)
    # Input values
    W_dg = InputValidate.get_value(df_aircraft, "MTOW", aircraft_idx)
    N_z = InputValidate.get_value(df_aircraft, "Load Factor Positive", aircraft_idx) * 1.5
    L = InputValidate.get_value(df_aircraft, "Fuselage Length", aircraft_idx)
    L_D = InputValidate.get_value(df_aircraft, "Fuselage Fineness Ratio", aircraft_idx) # REMEMBER! Not about L/D
    λ = InputValidate.get_value(df_aircraft,"Wing Taper Ratio",aircraft_idx) # Wing AR
    b_W = InputValidate.get_value(df_aircraft,"Wing Span",aircraft_idx) # Wing AR
    Λ = InputValidate.get_value(df_aircraft,"Wing Quarterchord Sweep",aircraft_idx) # Wing Quarterchord Sweep

    # Calculate wetted area from fuselage
    fuse = InputValidate.get_value(df_aircraft, "Fuselage Shape", aircraft_idx)
    ts = 0:0.01:1 # Discretisation
    S_f = wetted_area(fuse, ts) * u"m^2"

    # Assumptions
    K_door = 1.12  # Assume two side cargo doors
    K_Lg = 1.12  # Assume it is fuselage-mounted (if not, this is equal to 1.0)

    # Unit conversions
    W_dg = uconvert(u"lb", W_dg)
    L = uconvert(u"ft", L)
    S_f = uconvert(u"ft^2", S_f)
    b_W = uconvert(u"ft", b_W)

    # Estimate further values
    K_ws = ustrip(0.75 * ((1+2*λ)/(1+λ)) * b_W * tand(ustrip(Λ / L)))

    # Weight calculation
    W_fus = 0.3280 * K_door * K_Lg * (W_dg * N_z)^0.5 * L^0.25 * S_f^0.302 * (1.0 + K_ws)^0.04 * L_D^0.1

    return uconvert(u"kg", ustrip(W_fus) * u"lb")
end

function main_landing_gear_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    #TODO: CHECK! SEEMS RATHER LOW

    W_l = InputValidate.get_value(df_aircraft, "MLW", aircraft_idx)
    MTOW = InputValidate.get_value(df_aircraft, "MTOW", aircraft_idx)
    V_s = InputValidate.get_value(df_aircraft, "Landing Stall Speed", aircraft_idx)

    # Assumptions
    N_l = 1.5 * 3.0 # 1.5 x design load (3)
    K_mp = 1.0 # 1.126 for kneeling main gear 
    L_m = 3. * u"m" # Assume 3 m clearance to ground # TODO: best to approximate with better CG estimations to know tipback angle

    # Assume two structs to begin with
    N_mss = 2

    MTOW = uconvert(u"lb", MTOW)

    if MTOW < 200000.0*u"lb"
        N_mw = 4 # 2 x 2 
    else
        # Assume number of main wheels scale with 200,000 lb
        N_mw = ceil(ustrip(MTOW) / 200000.0) * 2 * 2

        # Roughly approximate how many struct needed (assume at most 6 wheels per struct)
        N_struct = ceil(N_mw / 6)

        # If more than number of structs, means each struct currently have more than 6 wheels,
        # take the new  N_struct
        if N_struct > N_mss
            N_mss = N_struct
        end
    end

    # Unit conversions
    W_l = uconvert(u"lb", W_l)
    V_s = uconvert(u"ft/s", V_s)
    L_m = uconvert(u"ft", L_m)

    W_mlg = 0.0106 * K_mp * W_l^0.888 * N_l^0.25 * L_m^0.4 * N_mw^0.321 * V_s^0.1 / N_mss^0.5

    return (uconvert(u"kg", ustrip(W_mlg) * u"lb"), N_mw)
end

function nose_landing_gear_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    #TODO: CHECK! SEEMS RATHER LOW

    W_l = InputValidate.get_value(df_aircraft, "MLW", aircraft_idx)

    # Assumptions
    N_l = 1.5 * 3.0 # 1.5 x design load (3)
    K_np = 1.0 # 1.126 for kneeling main gear 
    L_n = 3. * u"m" # Assume 3 m clearance to ground # TODO: best to approximate with better CG estimations to know tipback angle
    N_nw = 2 # Assume two tires

    # Unit conversions
    W_l = uconvert(u"lb", W_l)
    L_n = uconvert(u"ft", L_n)

    W_nlg = 0.032 * K_np * W_l^0.646 * N_l^0.2 * L_n^0.5 * N_nw^0.45

    return (uconvert(u"kg", ustrip(W_nlg) * u"lb"), N_nw)
end

function nacelle_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    K_ng = 1.017 # 1.0 for non-pylon mounted engine
    N_Lt = InputValidate.get_value(df_aircraft, "Engine Length Estimate", aircraft_idx)
    N_w = InputValidate.get_value(df_aircraft, "Engine Diameter Estimate", aircraft_idx)
    N_z = InputValidate.get_value(df_aircraft, "Load Factor Positive", aircraft_idx) * 1.5
    N_en = InputValidate.get_value(df_aircraft, "Number of Engines", aircraft_idx)

    # Calculate engine weight with content
    W_en = InputValidate.get_value(df_aircraft, "Engine Weight Estimate", aircraft_idx)
    engine_type = InputValidate.get_value(df_aircraft, "Engine Type", aircraft_idx)

    # Unit conversions
    N_Lt = uconvert(u"ft", N_Lt)
    N_w = uconvert(u"ft", N_w)
    W_en = uconvert(u"lb", W_en)
    
    if engine_type == "Jet"
        Kp = 1.0 # Propeller
        Ktr = 1.18 # Thrust reverser
    else
        Kp = 1.4 # Propeller
        Ktr = 1.0 # Thrust reverser
    end
    W_enc = 2.331 * Kp * Ktr * W_en^0.901

    # Get wetted area (using the cylinder model)
    engine = InputValidate.get_value(df_aircraft, "Engine Shape", aircraft_idx)
    ts = 0:0.01:1 # Discretisation
    S_n = wetted_area(engine[1],ts) * u"m^2" * N_en # Multiplied by the number of engine for total wetted area
    S_n = uconvert(u"ft^2", S_n)

    W_ncl = 0.6724 * K_ng * N_Lt^0.1 * N_w^0.294 * N_z^0.119 * W_enc^0.611 * N_en^0.984 * S_n^0.224

    return uconvert(u"kg", ustrip(W_ncl) * u"lb")
end

function engine_controls_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    N_en = InputValidate.get_value(df_aircraft, "Number of Engines", aircraft_idx)
    L_ec = InputValidate.get_value(df_aircraft, "Fuselage Length", aircraft_idx) # ASSUMPTION! Assume as long as the fuselage itself, roughly correct?

    # Unit conversions
    L_ec = uconvert(u"ft", L_ec)

    W_ec = 5 * N_en + 0.8 * ustrip(L_ec)

    return uconvert(u"kg", ustrip(W_ec) * u"lb")
end

function engine_starter_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    N_en = InputValidate.get_value(df_aircraft, "Number of Engines", aircraft_idx)
    W_en = InputValidate.get_value(df_aircraft, "Engine Weight Estimate", aircraft_idx)
    
    # Unit conversions
    W_en = uconvert(u"lb", W_en)

    W_estart = 49.19 * (N_en * ustrip(W_en) / 1000.0)^0.541

    return uconvert(u"kg", ustrip(W_estart) * u"lb")
end

function fuel_system_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    V_t = InputValidate.get_value(df_aircraft, "Fuel Tank Volume", aircraft_idx)
    
    # Unit conversions
    V_t = uconvert(u"gallon", V_t)

    # Assume 3 tanks... just an assumption!
    N_t = 3
    V_p = 0.5 * V_t
    V_i = 0.5 * V_t

    W_fs = 2.405 * V_t^0.606 * N_t^0.5 * (1 + V_p / V_t) / (1 + V_i / V_t)

    return uconvert(u"kg", ustrip(W_fs) * u"lb")
end

function flight_controls_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    # Assumptions
    N_f = 7 # Assume 7 flight functions
    N_m = 2 # Assume 2 mechanical controls

    # Get contorl surface estimate
    S_w = InputValidate.get_value(df_aircraft, "Wing Area", aircraft_idx)
    S_h = InputValidate.get_value(df_aircraft, "HT Area", aircraft_idx)
    S_t = InputValidate.get_value(df_aircraft, "VT Area", aircraft_idx)
    S_cs = 0.2 * (S_w+S_h+S_t) # Assume 20% of the area is for control surface
    S_cs = uconvert(u"ft^2", S_cs)

    # Get pitching inertia estimate
    L_ht = InputValidate.get_value(df_aircraft,"HT Tail Arm",aircraft_idx)  # Tail moment arm
    K_y = 0.3 * uconvert(u"ft", L_ht)
    MTOW = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx)  # MTOW
    MTOW = uconvert(u"lb", MTOW)
    I_y = MTOW * K_y^2

    W_fc = 145.9 * (N_f^0.554 * ustrip(S_cs)^0.2 * (ustrip(I_y) * 1e-6)^0.07) / (1 + N_m / N_f)

    return uconvert(u"kg", ustrip(W_fc) * u"lb")
end

function apu_installed_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    # TODO: Not going to implement this yet
    return 2.2 * W_apu
end

function instruments_weight(df_aircraft::DataFrame, aircraft_idx::Int,df_payload::DataFrame,payload_idx::Int)
    N_en = InputValidate.get_value(df_aircraft, "Number of Engines", aircraft_idx)
    L_f = InputValidate.get_value(df_aircraft, "Fuselage Length", aircraft_idx)
    B_w = InputValidate.get_value(df_aircraft, "Wing Span", aircraft_idx)

    # Get number of (flight) crew members
    N_c = InputValidate.get_value(df_payload, "Flight Crew", payload_idx)

    # Get engine type to decide on K_r and K_tp
    engine_type = InputValidate.get_value(df_aircraft, "Engine Type", aircraft_idx)

    # Default 1
    K_r = 1.0
    K_tp = 1.0

    if engine_type == "Turboprop"
        K_tp = 0.793
    elseif engine_type == "Propeller"
        K_r = 1.133
    end

    # Unit conversion
    B_w = uconvert(u"ft", B_w)
    L_f = uconvert(u"ft", L_f)

    W_inst = 4.509 * K_r * K_tp * N_c^0.541 * N_en * (L_f + B_w)^0.5

    return uconvert(u"kg", ustrip(W_inst) * u"lb")
end

function hydraulic_system_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    N_f = 7 # Assumption
    L_f = InputValidate.get_value(df_aircraft, "Fuselage Length", aircraft_idx)
    B_w = InputValidate.get_value(df_aircraft, "Wing Span", aircraft_idx)

    # Unit conversion
    B_w = uconvert(u"ft", B_w)
    L_f = uconvert(u"ft", L_f)

    W_hyd = 0.2673 * N_f * (L_f + B_w)^0.937

    return uconvert(u"kg", ustrip(W_hyd) * u"lb")
end

function electrical_system_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    #TODO: Not implemented yet
    return eee
end

function avionics_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    W_uav = 1400.0 # Assume uninstalled avionics is 1400 lbs, just upper end
    W_av = 1.73 * W_uav ^ 0.983

    return uconvert(u"kg", ustrip(W_av) * u"lb")
end

function furnishings_weight(df_aircraft::DataFrame, aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,payload_idx::Int)
    # Number of crew
    N_fc = InputValidate.get_value(df_payload, "Flight Crew", payload_idx)
    N_cc = InputValidate.get_value(df_payload, "Cabin Crew", payload_idx)
    N_pax = InputValidate.get_value(df_payload, "Passengers", payload_idx)
    N_p = N_pax + N_fc + N_cc

    # Seat Weight (in lbs)
    W_pax_seat = N_pax * 32.0 # Passenger
    W_fc_seat = N_fc * 60.0 # Flight Crew
    W_cc_seat = N_cc * 11.0 # Cabin Crew

    # Number of total crew
    N_c = N_fc + N_cc

    # Total cargo weight
    W_c = InputValidate.get_value(df_payload, "Cargo Weight", payload_idx) + InputValidate.get_value(df_payload, "Baggage Total Weight", payload_idx)
    
    # Calculate wetted area from fuselage
    fuse = InputValidate.get_value(df_aircraft, "Fuselage Shape", aircraft_idx)
    ts = 0:0.01:1 # Discretisation
    S_f = wetted_area(fuse, ts) * u"m^2"

    # Long range, short range assumptions
    range_idx = findfirst(==("Distance"), df_mission[:,1])
    max_range = maximum(skipmissing(collect(df_mission[range_idx,ncol(df_mission)-N_stages+1:ncol(df_mission)])))

    max_range = uconvert(u"km", max_range)

    short_haul_threshold = 500.0 * u"km"
    long_haul_threshold = 7000.0 * u"km"
    ultra_long_haul_threshold = 15000.0 * u"km"

    # This assumption neglects business jet furnishings
    if max_range < 500 * u"km"
        K_lav = 0.31
        K_buf = 1.02
    elseif max_range > ultra_long_haul_threshold
        K_lav = 1.11
        K_buf = 5.68
    elseif max_range > long_haul_threshold
        K_lav = 1.11
        K_buf = (5.68 - 1.02) / (ultra_long_haul_threshold - short_haul_threshold) * (max_range - short_haul_threshold) + 1.02
    else
        K_lav = (1.11 - 0.31) / (long_haul_threshold - short_haul_threshold) * (max_range - short_haul_threshold) + 0.31
        K_buf = (5.68 - 1.02) / (ultra_long_haul_threshold - short_haul_threshold) * (max_range - short_haul_threshold) + 1.02
    end

    # Unit conversion
    W_c = ustrip(uconvert(u"lb", W_c))
    S_f = ustrip(uconvert(u"ft^2", S_f))

    W_furn = 0.0577 * N_c^0.1 * W_c^0.393 * S_f^0.75 + (W_pax_seat + W_fc_seat + W_cc_seat) + K_lav * N_p^1.33 + K_buf * N_p^1.12

    return uconvert(u"kg", ustrip(W_furn) * u"lb")
end

function air_conditioning_weight(df_aircraft::DataFrame, aircraft_idx::Int,df_payload::DataFrame,payload_idx::Int)
    N_p = InputValidate.get_value(df_payload, "Passengers", payload_idx) + InputValidate.get_value(df_payload, "Flight Crew", payload_idx) + InputValidate.get_value(df_payload, "Cabin Crew", payload_idx)
    W_uav = 1400.0 # lbs assumption
    d = InputValidate.get_value(df_aircraft,"Diameter",aircraft_idx)
    l = InputValidate.get_value(df_aircraft,"Cabin Length",aircraft_idx) # Yes I am aware of the flight crew... No I am not going to include it
    V_pr = uconvert(u"ft^3", pi*(d/2)^2*l)

    W_airc = 62.36 * N_p^0.25 * (V_pr * 10^-3)^0.604 * W_uav^0.1

    return uconvert(u"kg", ustrip(W_airc) * u"lb")
end

function anti_icing_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    W_dg = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx) # Assume same as MTOW 
    W_dg = uconvert(u"lb",W_dg)
    W_ai = 0.002 * W_dg

    return uconvert(u"kg", W_ai)
end

function handling_gear_weight(df_aircraft::DataFrame, aircraft_idx::Int)
    regulations = InputValidate.get_value(df_aircraft,"Regulations",aircraft_idx)

    if occursin(r"(?i)Military",regulations)
        d = InputValidate.get_value(df_aircraft,"Diameter",aircraft_idx)
        l = InputValidate.get_value(df_aircraft,"Cabin Length",aircraft_idx)
        A_f = uconvert(u"ft^2", d*l)

        W_hg = ustrip(2.4 * A_f) * u"lb"
    else
        W_dg = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx) # Assume same as MTOW
        W_dg = uconvert(u"lb", W_dg)
        W_hg = 3.0 * (10^-4) * W_dg
    end

    return uconvert(u"kg", W_hg)
end

function fuel_calculation(df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    # Trapped fuel ratio (2%)
    trapped_fuel = InputValidate.get_value(df_aircraft,"Trapped Fuel Ratio",aircraft_idx)

    # Get old weight
    fuel_weight = InputValidate.get_value(df_aircraft,"Fuel Weight",aircraft_idx)
    MTOW = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx)

    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
    engine_type = InputValidate.get_value(df_aircraft,"Engine Type",aircraft_idx)

    # Estimation for weight remaining
    α = 1.0

    for col in ncol(df_mission)-N_stages+1:ncol(df_mission)
        stage = InputValidate.get_value(df_mission,"Stage",col)
        V = uconvert(u"m/s", InputValidate.get_value(df_mission,"Velocity",col))

        SFC_loiter = InputValidate.get_value(df_aircraft,"SFC_loiter",aircraft_idx)
        SFC_cruise = InputValidate.get_value(df_aircraft,"SFC_cruise",aircraft_idx)

        # SFC references
        if engine_type != "Jet"
            SFC_loiter = upreferred(SFC_loiter * V)
            SFC_cruise = upreferred(SFC_cruise * V)
        end

        if stage == "Takeoff"
            new_fraction = 0.97
        elseif stage == "Climb"
            new_fraction = 0.985
        elseif stage == "Loiter"
            LD = InputValidate.get_value(df_mission,"LD",col)
            E = uconvert(u"s", InputValidate.get_value(df_mission,"Duration",col))
            old_weight_fraction = InputValidate.get_value(df_mission,"Fuel Fraction",col)

            new_fraction = exp(-(E*SFC_loiter*g) / LD)
            df_mission = InputValidate.df_update_or_append(df=df_mission,label="Fuel Fraction",value=new_fraction,N_config=N_stages,col=col)
        elseif stage == "Cruise"
            LD = InputValidate.get_value(df_mission,"LD",col)
            R = uconvert(u"m", InputValidate.get_value(df_mission,"Distance",col))
            old_weight_fraction = InputValidate.get_value(df_mission,"Fuel Fraction",col)

            new_fraction = exp(-(R*SFC_cruise*g) / (LD * V))
            df_mission = InputValidate.df_update_or_append(df=df_mission,label="Fuel Fraction",value=new_fraction,N_config=N_stages,col=col)
        else
            new_fraction = 0.995
        end

        α *= new_fraction
        df_mission = InputValidate.df_update_or_append(df=df_mission,label="α",value=α,N_config=N_stages,col=col)
    end

    new_fuel_weight = MTOW*(1.0 - α)*trapped_fuel

    fuel_density = InputValidate.get_value(df_aircraft,"Fuel Density",aircraft_idx)
    
    fuel_tank_vol = (new_fuel_weight / fuel_density) * InputValidate.get_value(df_aircraft,"Fuel Volume Buffer",aircraft_idx)

    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Fuel Weight",value=new_fuel_weight,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Fuel Tank Volume",value=fuel_tank_vol,N_config=N_aircraft,col=aircraft_idx)

    return (df_aircraft, df_mission)
end

"""
    `weight_calculation` - A function which calculates a weight estimate

    return the weight estimation in dataframe
"""
function weight_calculation(;df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,N_payload::Int,payload_idx::Int)
    # Update fuel weight estimates
    (df_aircraft, df_mission) = fuel_calculation(df_aircraft,N_aircraft,aircraft_idx,df_mission,N_stages)
    
    # Aircraft wing weight
    W_wing = wing_weight(df_aircraft,aircraft_idx)
    
    # Horizontal Tail Weight
    W_HT = horizontal_tail_weight(df_aircraft,aircraft_idx)

    # Vertical Tail Weight
    W_VT = vertical_tail_weight(df_aircraft, aircraft_idx)

    # Fuselage Weight
    W_fus = fuselage_weight(df_aircraft, aircraft_idx)

    # Engine Weight
    W_eng = InputValidate.get_value(df_aircraft,"Engine Weight Estimate",aircraft_idx) * InputValidate.get_value(df_aircraft,"Number of Engines",aircraft_idx)
    W_ec = engine_controls_weight(df_aircraft, aircraft_idx)
    W_estart = engine_starter_weight(df_aircraft, aircraft_idx)
    W_ncl = nacelle_weight(df_aircraft, aircraft_idx) # Nacelle Weight

    # Landing Gear Weight
    (W_mlg, N_mw) = main_landing_gear_weight(df_aircraft, aircraft_idx)
    (W_nlg, N_nw) = nose_landing_gear_weight(df_aircraft, aircraft_idx)

    # Numebr of landing gear tires
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Number of Tires",value=N_mw+N_nw,N_config=N_aircraft,col=aircraft_idx)

    # Fuel System
    W_fs = fuel_system_weight(df_aircraft, aircraft_idx)

    # Flight Control System
    W_fc = flight_controls_weight(df_aircraft, aircraft_idx)

    # Instruments
    W_inst = instruments_weight(df_aircraft, aircraft_idx,df_payload,payload_idx)

    # Hydraulics
    W_hyd = hydraulic_system_weight(df_aircraft, aircraft_idx)

    # Avionics
    W_av = avionics_weight(df_aircraft, aircraft_idx)

    # Furnishing
    W_furn = furnishings_weight(df_aircraft, aircraft_idx,df_mission,N_stages,df_payload,payload_idx)

    # Air Conditioning
    W_airc = air_conditioning_weight(df_aircraft, aircraft_idx,df_payload,payload_idx)

    # Anti-ice
    W_ai = anti_icing_weight(df_aircraft, aircraft_idx)

    # Handling gear
    W_hg = handling_gear_weight(df_aircraft, aircraft_idx)

    W_total = W_wing + W_HT + W_VT + W_fus + W_eng + W_mlg + W_nlg + W_ncl + W_ec + W_estart + W_fs + W_fc + W_inst + W_hyd + W_av + W_furn + W_airc + W_ai + W_hg
    print("Current: ")
    print(W_total)
    print("\nPredicted: ")
    print(InputValidate.get_value(df_aircraft,"Empty Weight",aircraft_idx))
    print("\n")

    # Obtain fuel and paylaod weight to calculate the MTOW
    W_fuel = InputValidate.get_value(df_aircraft,"Fuel Weight",aircraft_idx)
    W_pl = InputValidate.get_value(df_payload,"Payload Weight",payload_idx)

    MTOW = W_total + W_fuel + W_pl

    print("Current MTOW: ")
    print(MTOW)
    print("\nPredicted: ")
    print(InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx))
    print("\n")

    # Update the weight values with new ones
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Empty Weight",value=W_total,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="MTOW",value=W_total,N_config=N_aircraft,col=aircraft_idx)

    return df_aircraft
end

end