using Plots
using Unitful

include("constraintdiagram.jl")

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

sample_engine = ConstraintDiagram.Engine(
    type = "Jet",
    λ_bpr = 10,
    η_prop = 0.8,
    prop_type = "Constant Speed",
    prop_disc_load = 100
)

sample_aircraft = ConstraintDiagram.Aircraft(
    engine = sample_engine,
    CL_max = 1.5, # CL_max = Maximum coefficient of lift
    CD0 = 0.02, # CD0 = Zero-lift drag
    AR  = 10,   # AR = Aspect ratio
    e   = 0.8,    # e = Oswald Efficiency
    N_e = 4
)

sample_point = ConstraintDiagram.PointPerformance(
    aircraft   = sample_aircraft,
    α          = 1, # α = current weight / total weight ratio
    β          = 1, # β = thrust / sea-level thrust ratio
    h          = 25000, # h = height of aircraft (default h_unit = ft)
    V∞         = 200 # V∞ = Freestream velocity (default V∞_unit = knots)
)

sample_climb = ConstraintDiagram.PointPerformance(
    aircraft   = sample_aircraft,
    α          = 1, # α = current weight / total weight ratio
    β          = 1, # β = thrust / sea-level thrust ratio
    h          = 500, # h = height of aircraft (default h_unit = ft)
    V∞         = 170, # V∞ = Freestream velocity (default V∞_unit = knots)
    G          = 2.4,
    N_a        = 3
)

sample_takeoff = ConstraintDiagram.TakeOffPerformance(
    aircraft   = sample_aircraft,
    TODA       = 2119, # Takeoff distance
    TODA_unit  = "m",
    α          = 1, # α = current weight / total weight ratio
    β          = 1, # β = thrust / sea-level thrust ratio
    h          = 0, # h = height of aircraft (default h_unit = ft)
    surface    = "Concrete" # surface = Surface for takeoff (concrete or grass for now)
)

sample_landing = ConstraintDiagram.LandingPerformance(
    aircraft   = sample_aircraft,
    LDA        = 2119, # Landing distance
    LDA_unit  = "m",
    α          = 0.87, # α = current weight / total weight ratio
    h          = 0, # h = height of aircraft (default h_unit = ft)
    K_r        = 1
)

(WS_to, TW_to) = ConstraintDiagram.bfl_perf(sample_takeoff, 100, 5000, 1000)

WS_land = ConstraintDiagram.landing_perf(sample_landing)

(WS_c, TW_c) = ConstraintDiagram.point_perf(sample_point, 100, 5000, 1000)

(WS_climb, TW_climb) = ConstraintDiagram.point_perf(sample_climb, 100, 5000, 1000)


## Plotting
gr()
p = plot(WS_to, TW_to, label = "BFL")
plot!(p, WS_c, TW_c, label = "Cruise")
plot!(p, WS_climb, TW_climb, label = "Climb")
vline!([WS_land], label = "Landing")
title!(p, "Constraint Diagram")
xlabel!(p, join(["Wing Loading (W/S) ", u"kg / m / s ^ 2"]))

if sample_engine.type == "Jet"
    ylabel!(p, "Thrust to Weight Ratio (T/W)")
    ylims!(p, 0, 1)
else
    ylabel!(p, join(["Power to Weight Ratio (P/W) ", u"m / s"]))
end

