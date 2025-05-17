module AircraftAero

# Initialise packages used
using DataFrames
using Unitful
using AeroFuse
using Roots
using LinearAlgebra

# Include the InputValidate module
include("inputvalidate.jl")

"""
    `make_case` - A function which runs VLM

    (From AeroFuse Examples)

    return the VLM analysis
"""
function make_case(α, wing_mesh, refs)
    aircraft = ComponentVector(wing = make_horseshoes(wing_mesh))

    # Freestream conditions
    fs = Freestream(alpha = α)  # Design variable: Angle of attack

    # Solve system
    return VortexLatticeSystem(aircraft, fs, refs, true)
end

"""
    `get_forces` - A function which gets the forces from an analysis

    (From AeroFuse Examples)

    return the forces
"""
function get_forces(system, wing_mesh)
    # Evaluate aerodynamic coefficients
    CDi, CY, CL, Cl, Cm, Cn = nearfield(system)
    # CDi, _, _ = farfield(system)

    # Calculate equivalent flat-plate skin-friction drag
    # CDv = parasitic_drag_coefficient(wing_mesh, 1.0, system.reference)

    # Calculate local-dissipation/local-friction drag
    CVs = norm.(surface_velocities(system)).wing
    CDv = parasitic_drag_coefficient(wing_mesh, system.reference, 0.8, CVs)

    return (CDi = CDi, CDv = CDv, CD = CDi + CDv, CL = CL)
end

function initiate_wing_mesh(df_aircraft,aircraft_idx, wing_type,n)
    # Get data
    taper_ratio = InputValidate.get_value(df_aircraft,"$wing_type Taper Ratio",aircraft_idx)
    b = InputValidate.get_value(df_aircraft,"$wing_type Span",aircraft_idx)
    S = InputValidate.get_value(df_aircraft,"$wing_type Area",aircraft_idx)
    aerofoil = InputValidate.get_value(df_aircraft,"$wing_type Airfoil",aircraft_idx)
    sweep = InputValidate.get_value(df_aircraft,"$wing_type Quarterchord Sweep",aircraft_idx)
    root_twist = InputValidate.get_value(df_aircraft,"$wing_type Setting Angle",aircraft_idx)
    tip_twist = InputValidate.get_value(df_aircraft,"$wing_type Twist Angle",aircraft_idx)
    dihedral = InputValidate.get_value(df_aircraft,"$wing_type Dihedral",aircraft_idx)

    # Get chord
    c_root = 2*S/(b*(1+taper_ratio))
    c_tip = taper_ratio * c_root
    chords = LinRange(c_root, c_tip, n)

    # Get halfsapn
    bs = fill(b/2/(n-1), n-1)

    # Get twist, dihedral and sweep
    twists = LinRange(root_twist, tip_twist, n)
    dihedrals = fill(dihedral, n-1)
    sweeps = fill(sweep, n-1)

    if occursin(r"(?i)^NACA\d{4}$", aerofoil)
        x1 = aerofoil[5]
        x2 = aerofoil[6]
        x3 = aerofoil[7]
        x4 = aerofoil[8]
        foil = fill(naca4(x1,x2,x3,x4))
    else
        foil = AeroFuse.AircraftGeometry.read_foil("./airfoil_database/$aerofoil.dat") 
    end

    # Create wing
    wing = Wing(
        foils     = fill(foil, n),
        chords    = ustrip.(chords), # Design variables
        twists    = ustrip.(twists),
        spans     = ustrip.(bs), # Normalizing halfspan to 1
        dihedrals = ustrip.(dihedrals),
        sweeps    = ustrip.(sweeps), # Quarter-chord sweep
        w_sweep   = 0.25,
        symmetry  = true
    )

    # Meshing
    wing_mesh = WingMesh(
        wing, fill(2, n - 1), 6, 
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
    MTOW = InputValidate.get_value(df_aircraft,"MTOW",aircraft_idx)
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
    n_vars = 64 # Number of discretisation

    # HT and VT sizing
    

    # Wing calculation
    (wing, wing_mesh) = initiate_wing_mesh(df_aircraft,aircraft_idx,"Wing",n_vars)
    Sw = projected_area(wing_mesh) # Reference area

    # HT calculation

    row_idx = findfirst(==("Stage"),df_mission[:,1])
    cruise_condition = findall(==("Cruise"),skipmissing(collect(df_mission[row_idx, :])))

    for col in cruise_condition
        alpha = InputValidate.get_value(df_mission,"α",col)
        ρ = InputValidate.get_value(df_mission,"ρ",col)
        V = uconvert(u"m/s",InputValidate.get_value(df_mission,"Velocity",col))

        # Get target CL at cruise
        CL_cruise = ustrip((MTOW*alpha*g) / (0.5*ρ*V^2*Sw))  # Target lift coefficient
        print("CL_cruise: \n")
        print(CL_cruise)
        print("\n")

        refs = References(
            speed     = ustrip(V),
            area      = Sw,
            density   = ustrip(ρ),
            span      = span(wing_mesh),
            chord     = mean_aerodynamic_chord(wing_mesh),
            location  = mean_aerodynamic_center(wing_mesh)
        )

        # Find angle of attack which matches target CL
        α0 = find_zero(4.0, Roots.Order0()) do α
            sys = make_case(α, wing_mesh, refs)
            CL_cruise - get_forces(sys, wing_mesh).CL
        end

        print("α0:\n")
        print(α0)
        print("\n")

        sys = make_case(α0, wing_mesh, refs)
        init = get_forces(sys, wing_mesh)
        print_coefficients(sys)
        print(init)
    end

    return df_aircraft
end

end