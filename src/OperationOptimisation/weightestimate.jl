module WeightEst

# Initialise packages used
using DataFrames
using Unitful
using AeroFuse

# Include the InputValidate module
include("inputvalidate.jl")

# Include AircraftAero module
include("aircraftaero.jl")

function wing_weight(df_aircraft,aircraft_idx)
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

function horizontal_tail_weight(df_aircraft, aircraft_idx)
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

function vertical_tail_weight(df_aircraft, aircraft_idx)
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


"""
    `weight_calculation` - A function which calculates a weight estimate

    return the weight estimation in dataframe
"""
function weight_calculation(;df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,N_payload::Int,payload_idx::Int)
    # Aircraft wing weight
    W_wing = wing_weight(df_aircraft,aircraft_idx)
    
    # Horizontal Tail Weight
    W_HT = horizontal_tail_weight(df_aircraft,aircraft_idx)

    # Vertical Tail Weight
    W_VT = vertical_tail_weight(df_aircraft, aircraft_idx)

    # Fuselage Weight
    W_fus = fuselage_weight(df_aircraft, aircraft_idx)

    # Engine Weight
    W_eng = InputValidate.get_value(df_aircraft,"Engine Weight Estimate",aircraft_idx)

    W_total = W_wing + W_HT + W_VT + W_fus + W_eng
    print("Current: ")
    print(W_total)
    print("\nPredicted: ")
    print(InputValidate.get_value(df_aircraft,"Empty Weight",aircraft_idx))
    print("\n")


    return df_aircraft
end

end