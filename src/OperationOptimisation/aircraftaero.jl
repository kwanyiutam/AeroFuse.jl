module AircraftAero

# Initialise packages used
using DataFrames
using Unitful
using AeroFuse
using Roots
using LinearAlgebra
using Plots

# Include the InputValidate module
include("inputvalidate.jl")


#### Additional functions for engine (nacelles)
parasitic_drag_coefficient_eng(eng :: HyperEllipseFuselage, refs :: References, x_tr :: Real, ts = 0:0.01:1) = wetted_area_drag_coefficient(eng, x_tr, refs.density, refs.speed, mach_number(refs), refs.viscosity, refs.area, ts)

form_factor_engine(f) = 1 + 0.35/f

form_factor(eng :: HyperEllipseFuselage) = form_factor_engine(eng.length / eng.radius)

function wetted_area_drag_coefficient(eng :: HyperEllipseFuselage, x_tr, ρ, V, M, μ, S_ref, ts = 0:0.01:1)
    # Fuselage quantities
    L = eng.length
    S_wet = wetted_area(eng, ts) 
    Kf = form_factor(eng)
    fM = (1 - 0.08M^1.45)

    # Parasitic drag coefficient
    CDp = AeroFuse.parasitic_drag_coefficient(L, x_tr, ρ, V, M, μ, S_ref, S_wet, Kf, fM)

    return CDp
end

"""
    `make_case` - A function which runs VLM

    (From AeroFuse Examples)

    return the VLM analysis
"""
function make_case(α, wing_mesh, HT_mesh, VT_mesh, refs)
    aircraft = ComponentVector(
        wing = make_vortex_rings(wing_mesh),
        htail = make_vortex_rings(HT_mesh),
        vtail = make_vortex_rings(VT_mesh)
    )

    fs = Freestream(alpha = α)  # Design variable: Angle of attack

    # Solve system
    return VortexLatticeSystem(aircraft, fs, refs, true)
end

"""
    `get_forces` - A function which gets the forces from an analysis

    (From AeroFuse Examples)

    return the forces
"""
function get_forces(system, wing, HT, VT, fuse, eng_save, CD_upsweep)
    # Evaluate aerodynamic coefficients
    CDi, CY, CL, Cl, Cm, Cn = nearfield(system)
    # CDi, _, _ = farfield(system)

    # Calculate equivalent flat-plate skin-friction drag
    # CDv = parasitic_drag_coefficient(wing_mesh, 1.0, system.reference)

    # Wing viscous drag
    CVs = norm.(surface_velocities(system)).wing
    CDv = parasitic_drag_coefficient(wing, system.reference, 0.1, CVs)

    # HT viscous drag
    CVs_HT = norm.(surface_velocities(system)).htail
    CDv_HT = parasitic_drag_coefficient(HT, system.reference, 0.1, CVs_HT)

    # VT viscous drag
    CVs_VT = norm.(surface_velocities(system)).vtail
    CDv_VT = parasitic_drag_coefficient(VT, system.reference, 0.1, CVs_VT)

    # Fuselage viscous drag
    CDv_F = parasitic_drag_coefficient(fuse, system.reference, 0.05)

    CDv_E = 0
    for i in eachindex(eng_save)
        CDv_E = CDv_E + parasitic_drag_coefficient_eng(eng_save[i], system.reference, 0.001)
    end

    CDv = CDv + CDv_HT + CDv_VT + CDv_F + CDv_E + CD_upsweep

    return (CDi = CDi, CDv = CDv, CDv_HT = CDv_HT, CDv_VT = CDv_VT, CDv_F = CDv_F, CD = CDi + CDv, CL = CL, L_D = CL / (CDi + CDv))
end

function initiate_wing_mesh(df_aircraft,aircraft_idx,wing_type,n_span,n_chord)
    # Define
    if wing_type == "Wing"
        angle = 0.
        axis = [0., 1., 0.]
        symmetry = true
    elseif wing_type == "HT"
        angle = 0.
        axis = [0., 1., 0.]
        symmetry = true
    elseif wing_type == "VT"
        angle = 90.
        axis = [1., 0., 0.]
        symmetry = false
    else
        @warn "wing_type $wing_type cannot be found, assuming it is the same orientation as wing..."
        angle = 0.
        axis = [0., 0., 0.]
        symmetry = true
    end

    # Get data
    positions = InputValidate.get_value(df_aircraft,"$wing_type Position",aircraft_idx)
    taper_ratio = InputValidate.get_value(df_aircraft,"$wing_type Taper Ratio",aircraft_idx)
    b = uconvert(u"m", InputValidate.get_value(df_aircraft,"$wing_type Span",aircraft_idx))
    S = uconvert(u"m^2", InputValidate.get_value(df_aircraft,"$wing_type Area",aircraft_idx))
    aerofoil = InputValidate.get_value(df_aircraft,"$wing_type Airfoil",aircraft_idx)
    sweep = InputValidate.get_value(df_aircraft,"$wing_type Quarterchord Sweep",aircraft_idx)
    root_twist = InputValidate.get_value(df_aircraft,"$wing_type Setting Angle",aircraft_idx)
    tip_twist = InputValidate.get_value(df_aircraft,"$wing_type Twist Angle",aircraft_idx)
    dihedral = InputValidate.get_value(df_aircraft,"$wing_type Dihedral",aircraft_idx)

    # Number of spanwsie changes specified
    n_change = 2

    # Get chord
    c_root = 2*S/(b*(1+taper_ratio))
    c_tip = taper_ratio * c_root
    chords = LinRange(c_root, c_tip, n_change)

    # Get halfsapn
    bs = fill(b/2/(n_change-1), n_change-1)

    # Get twist, dihedral and sweep
    twists = LinRange(root_twist, tip_twist, n_change)
    dihedrals = fill(dihedral, n_change-1)
    sweeps = fill(sweep, n_change-1)

    if occursin(r"(?i)^NACA\d{4}$", aerofoil)
        x1 = parse(Int, aerofoil[5])
        x2 = parse(Int, aerofoil[6])
        x3 = parse(Int, aerofoil[7])
        x4 = parse(Int, aerofoil[8])
        foil = naca4(x1,x2,x3,x4)
    else
        foil = AeroFuse.AircraftGeometry.read_foil("./airfoil_database/$aerofoil.dat") 
    end

    # Create wing
    wing = Wing(
        foils     = fill(foil, n_change),
        chords    = ustrip.(chords), # Design variables
        twists    = ustrip.(twists),
        spans     = ustrip.(bs),
        dihedrals = ustrip.(dihedrals),
        sweeps    = ustrip.(sweeps), # Quarter-chord sweep
        position  = ustrip.(positions),
        angle     = angle,
        axis      = axis,
        w_sweep   = 0.25,
        symmetry  = symmetry
    )

    # Meshing
    wing_mesh = WingMesh(
        wing, fill(n_span,n_change-1), n_chord, 
        span_spacing = Cosine(),
    )

    return (wing, wing_mesh)
end

"""
    `run_aero_analysis` - A function which runs the aerodynamics analysis for OperationOptimisation

    return the original dataframe
"""
function run_aero_analysis(;df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int)
    # Key variables
    MTOW = uconvert(u"kg", InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx))
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
    n_span = 20 # Number of spanwise discretisation
    n_chord = 8 # Number of chordwise discretisation

    # Obtain wing position
    fuselage_length = uconvert(u"m", InputValidate.get_value(df_aircraft,"Fuselage Length",aircraft_idx))
    diameter = uconvert(u"m", InputValidate.get_value(df_aircraft,"Diameter",aircraft_idx))
    wing_x_pos = InputValidate.get_value(df_aircraft,"Wing MAC Location (Relative to Fuselage Length)",aircraft_idx) * fuselage_length
    wing_z_pos = InputValidate.get_value(df_aircraft,"Wing Root Chord Location (Relative to Fuselage Centerline, Normalised by Radius)",aircraft_idx) * diameter / 2
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Wing Position",value=(wing_x_pos,0.0*u"m",wing_z_pos),N_config=N_aircraft,col=aircraft_idx)

    # Wing calculation
    (wing, wing_mesh) = initiate_wing_mesh(df_aircraft,aircraft_idx,"Wing",n_span,n_chord)
    Sw = projected_area(wing_mesh) # Reference area    

    ###HT and VT calculation
    # Get the key values (dimensions and volumetric coefficient)
    c_mean = mean_aerodynamic_chord(wing_mesh) *u"m"
    bref = span(wing_mesh) * u"m"
    HT_Vbar = InputValidate.get_value(df_aircraft,"VH",aircraft_idx)
    VT_Vbar = InputValidate.get_value(df_aircraft,"VV",aircraft_idx)

    # Get the dimensions
    HT_x_pos = InputValidate.get_value(df_aircraft,"HT MAC Location (Relative to Fuselage Length)",aircraft_idx) * fuselage_length
    HT_AR = InputValidate.get_value(df_aircraft,"HT AR",aircraft_idx)
    VT_x_pos = InputValidate.get_value(df_aircraft,"VT MAC Location (Relative to Fuselage Length)",aircraft_idx) * fuselage_length
    VT_AR = InputValidate.get_value(df_aircraft,"VT AR",aircraft_idx)

    # Calculate HT and VT properties
    S_HT = HT_Vbar * Sw*u"m^2" * c_mean / (HT_x_pos-wing_x_pos)
    b_HT = sqrt(HT_AR * S_HT)
    S_VT = VT_Vbar * Sw*u"m^2" * bref / (VT_x_pos-wing_x_pos)
    b_VT = sqrt(VT_AR * S_VT)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="HT Area",value=S_HT,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="HT Span",value=b_HT,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="VT Area",value=S_VT,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="VT Span",value=b_VT,N_config=N_aircraft,col=aircraft_idx)

    # Calcualte vertical position (based on the assumption that the tail follows the afterbody upsweep)
    HT_z_pos = InputValidate.get_value(df_aircraft,"HT Root Chord Location (Relative to Span of VT)",aircraft_idx) * b_VT
    tail_angle = InputValidate.get_value(df_aircraft,"Afterbody Upsweep Angle",aircraft_idx)
    nose_x_end = InputValidate.get_value(df_aircraft,"Nose Length",aircraft_idx)
    tail_x_start = nose_x_end + InputValidate.get_value(df_aircraft,"Cabin Length",aircraft_idx)
    VT_z_pos = (VT_x_pos-tail_x_start)*tand(tail_angle) - (diameter / 2)
    HT_z_pos = HT_z_pos + VT_z_pos

    # Save data into dataframe
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="HT Position",value=(HT_x_pos,0.0*u"m",HT_z_pos),N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="VT Position",value=(VT_x_pos,0.0*u"m",VT_z_pos),N_config=N_aircraft,col=aircraft_idx)

    # HT and VT calculation
    (HT, HT_mesh) = initiate_wing_mesh(df_aircraft,aircraft_idx,"HT",n_span,n_chord)
    (VT, VT_mesh) = initiate_wing_mesh(df_aircraft,aircraft_idx,"VT",n_span,n_chord)

    # Initiate Fuselage
    nose_end = nose_x_end / fuselage_length
    rear_start = tail_x_start / fuselage_length
    fuse = HyperEllipseFuselage(
        radius = ustrip(diameter) / 2.0,
        length = ustrip(fuselage_length),
        x_a = ustrip(nose_end),
        x_b = ustrip(rear_start),
        d_rear = ustrip(VT_z_pos),
        c_nose = 1.5,
        c_rear = 1.5,
    )

    # Initiate Nacelles
    N_engines = InputValidate.get_value(df_aircraft,"Number of Engines",aircraft_idx)
    engine_diameter = uconvert(u"m", InputValidate.get_value(df_aircraft,"Engine Diameter Estimate",aircraft_idx))
    engine_length = uconvert(u"m", InputValidate.get_value(df_aircraft,"Engine Length Estimate",aircraft_idx))
    engine_mount = InputValidate.get_value(df_aircraft,"Engine Mounting Position",aircraft_idx)
    engine_sep = InputValidate.get_value(df_aircraft,"Engine Position from Fuselage Normalised by Engine Diameter",aircraft_idx)
    engine_z_pos = wing_z_pos - engine_diameter

    if engine_mount == "Wing"
        engine_x_pos = wing_x_pos
    elseif engine_mount == "Tail"
        engine_x_pos = fuselage_length
    else
        throw(ArgumentError("Invalid Argument, cannot find the engine mount position of $engine_mount"))
    end

    eng_save = []
    N_engine_half = floor(N_engines/2)
  
    for i in 1:N_engines
        engine_y_pos = -diameter/2.0 -((N_engine_half-i+1)*engine_sep+0.5+N_engine_half-i)*engine_diameter # Fuselage and engine separated by two diameter length

        if i > N_engine_half
            engine_y_pos = engine_y_pos + diameter + engine_sep*engine_diameter
        end

        if (i == N_engines) && (N_engines % 2 == 1)
            engine_x_pos = fuselage_length
            engine_y_pos = 0.0
            engine_z_pos = VT_z_pos
        end

        position = [engine_x_pos,engine_y_pos,engine_z_pos]

        eng = HyperEllipseFuselage(
            radius = ustrip(engine_diameter) / 2.0,
            length = ustrip(engine_length),
            x_a = 0.1,
            x_b = 0.8,
            c_nose = 2,
            c_rear = 2,
            position = ustrip.(position),
        )

        push!(eng_save, eng)
    end

    #gr()

    ## Coordinates
    #Plots.plot(
    #    aspect_ratio = 1,
    #    camera = (30, 30),
    #    zlim = span(wing) .* (-0.5, 0.5),
    #    size = (800, 600)
    #)
    #Plots.plot!(wing_mesh, label = "Wing")
    #Plots.plot!(HT_mesh, label = "HT")
    #Plots.plot!(VT_mesh, label = "VT")
    #Plots.plot!(fuse, label = "Fuselage")

    #for i in 1:N_engines
    #    Plots.plot!(eng_save[i], label = "Engine $i")
    #end

    #savefig("SamplePlane.png") 

    # Add upsweep related drag
    CD_upsweep = ustrip(3.83*(pi*diameter^2/(4*Sw))*deg2rad(tail_angle)^2.5)

    # Get all conditions for cruise/loiter
    row_idx = findfirst(==("Stage"),df_mission[:,1])
    cruise_condition = findall(==("Cruise"),skipmissing(collect(df_mission[row_idx, :])))
    loiter_condition = findall(==("Loiter"),skipmissing(collect(df_mission[row_idx, :])))
    cruiseloiter_condition = vcat(cruise_condition,loiter_condition)

    for col in cruiseloiter_condition
        alpha = InputValidate.get_value(df_mission,"α",col)
        ρ = InputValidate.get_value(df_mission,"ρ",col)
        V = uconvert(u"m/s",InputValidate.get_value(df_mission,"Velocity",col))

        # Get target CL at cruise
        CL_cruise = ustrip((MTOW*alpha*g) / (0.5*ρ*V^2*Sw))  # Target lift coefficient

        refs = References(
            speed     = ustrip(V),
            area      = Sw,
            density   = ustrip(ρ),
            span      = span(wing_mesh),
            chord     = mean_aerodynamic_chord(wing_mesh),
            location  = mean_aerodynamic_center(wing_mesh)
        )

        # Find angle of attack which matches target CL
        α0 = find_zero(2.0, Roots.Order0()) do α
            sys = make_case(α, wing_mesh, HT_mesh, VT_mesh, refs)
            CL_cruise - get_forces(sys, wing_mesh, HT_mesh, VT_mesh, fuse, eng_save, CD_upsweep).CL
        end

        sys = make_case(α0, wing_mesh, HT_mesh, VT_mesh, refs)
        init = get_forces(sys, wing_mesh, HT_mesh, VT_mesh, fuse, eng_save, CD_upsweep)

        df_mission = InputValidate.df_update_or_append(df=df_mission,label="LD",value=init.L_D,N_config=N_stages,col=col)
        df_mission = InputValidate.df_update_or_append(df=df_mission,label="CD",value=init.CD,N_config=N_stages,col=col)
        df_mission = InputValidate.df_update_or_append(df=df_mission,label="CD0_Estimate",value=init.CDv,N_config=N_stages,col=col)
    end

    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Wing Mesh",value=wing_mesh,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="HT Mesh",value=HT_mesh,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="VT Mesh",value=VT_mesh,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Fuselage Shape",value=fuse,N_config=N_aircraft,col=aircraft_idx)
    df_aircraft = InputValidate.df_update_or_append(df=df_aircraft,label="Engine Shape",value=eng_save,N_config=N_aircraft,col=aircraft_idx)

    return (df_aircraft, df_mission)
end

end