# This is the main constraint diagram file, consisting of
# all the required modules and packages
# some functions (e.g. aircraftdefine) are called instead of defined here

module ConstraintDiagram

# Initialise packages used
using ISAData
using Unitful
using DataFrames
using Interpolations

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the ParseData module
include("parsedata.jl")

# Include the AircraftDefine functions
include("aircraftdefine.jl")

# Include the MissionDefine functions
include("missiondefine.jl")

# INclude the weight cost functions
include("weightcostcalc.jl")

##### Point Performance Definition
abstract type AbstractPointPerformance end

struct PointPerformance{T} <: AbstractPointPerformance
    aircraft    :: Aircraft{T}
    N_a         :: T # N_a = Number of available engines
    α           :: Array{T} # α = current weight / total weight ratio
    β           :: T # β = thrust / sea-level thrust ratio
    h           :: T # h = height of aircraft (default h_unit = ft)
    h_unit      :: String 
    V∞          :: Array{T} # V∞ = Freestream velocity (default V∞_unit = knots)
    V∞_unit     :: String
    G           :: T # G = climb gradient (%)
    a_max       :: T # a_max = maximum acceleration (default a_max_unit = m/s^2)
    a_max_unit  :: String 
    n           :: T # n = Load factor
end

function PointPerformance(; aircraft :: Aircraft, N_a :: Int = 2, α, β, h, h_unit :: String = "ft", V∞, V∞_unit :: String = "knots", G = 0, a_max = 0, a_max_unit :: String = "m/s^2", n = 1)
    T = promote_type(eltype(aircraft.N_e), eltype(aircraft.CD0), eltype(aircraft.e), eltype(aircraft.engine.η_prop), eltype(aircraft.engine.λ_bpr), eltype(N_a), eltype(α), eltype(β), eltype(h), eltype(V∞), eltype(G), eltype(a_max), eltype(n))
    @assert N_a > 0 && aircraft.N_e >= N_a "Number of available engines must be bigger than 1 and less than or equal to number of engines."
    @assert h >= 0 && all(>=(0), any(isnan, V∞)) "h and V∞ inputs must be non-negative"
    @assert all(>=(0), any(isnan, α)) && all(<=(1), any(isnan, α)) && β >= 0 && β <= 1 "α and β must be between 0 and 1"
    
    return PointPerformance{T}(aircraft, N_a, α, β, h, h_unit, V∞, V∞_unit, G, a_max, a_max_unit, n)
end

## Constraint Diagram Build-up
#=======================#

# Point performance calculation: 
# WS_min = Minimum Wing Loading considered
# WS_max = Maximum Wing Loading considered
# WS_N = Number of Wing Loading points
function point_perf(
    point_data :: PointPerformance,
    WS_min = 1,
    WS_max = 10000,
    WS_N = 1000,
    )

    ### Convert units
    # Convert to unit form, especially for knots
    h_unit = AeroUnits.convert_to_unit(point_data.h_unit)
    V∞_unit = AeroUnits.convert_to_unit(point_data.V∞_unit)
    a_max_unit = AeroUnits.convert_to_unit(point_data.a_max_unit)

    # Converting to the correct units
    h = uconvert(u"m", point_data.h * h_unit)
    V∞ = uconvert.(u"m/s", point_data.V∞ .* V∞_unit)
    a_max = uconvert(u"m/s^2", point_data.a_max * a_max_unit)

    η_prop = point_data.aircraft.engine.η_prop
    α = point_data.α

    ### Determine whether the input is power or thrust 
    if point_data.aircraft.engine.type == "Jet"
        prop_constant = 1
    else
        prop_constant = V∞ / η_prop
    end

    ### Determine density
    (ρ,_,_,_) = ISAdata(h)

    ### Terms in the point-performance
    WS = LinRange(WS_min*u"kg/m/s^2", WS_max*u"kg/m/s^2", WS_N)
    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
    C_climb = sin(atan(point_data.G / 100)) # Related to climb gradient
    C_speed = a_max / g # Related to change in speed
    C_zerol = (0.5 .* ρ .* (V∞ .^ 2) .* point_data.aircraft.CD0) ./ (α .* WS)
    C_induc = (α .* point_data.n ^2 .* WS) ./ (0.5 .* ρ .* (V∞ .^ 2) .* pi .* point_data.aircraft.AR .* point_data.aircraft.e)

    TW = (α ./ point_data.β) .* (point_data.aircraft.N_e ./ point_data.N_a) .* prop_constant .* (C_climb .+ C_speed .+ C_zerol .+ C_induc)

    return (WS, TW)
end

##### Takeoff Performance Definition
abstract type AbstractTakeOffPerformance end

struct TakeOffPerformance{T} <: AbstractTakeOffPerformance
    aircraft    :: Aircraft{T}
    TODA        :: T # TODA = Take Off Distance Available (default TODA_unit = ft)
    TODA_unit   :: String
    α           :: Array{T} # α = current weight / total weight ratio
    β           :: T # β = thrust / sea-level thrust ratio
    h           :: T # h = height of aircraft (default h_unit = ft)
    h_unit      :: String 
    γ_2         :: T # γ_2 = Second segment climb
    γ_2_min     :: T # γ_2 = Second segment climb (minimum, based on engine)
    h_to        :: T # h_to = clearance height (35 ft for FAR and 50 ft for military) (default h_to_unit = ft)
    h_to_unit   :: String
    V_2_coef    :: T # V_2_coef = V_2 speed as a function of stall speed, usually 1.2
    μ           :: T # μ = Coefficient of friction of ground
end

function TakeOffPerformance(; aircraft :: Aircraft, TODA, TODA_unit :: String = "ft", α, β, h = 0, h_unit :: String = "ft", γ_2_min = 0.24, γ_2 = γ_2_min, surface :: String = "Concrete", h_to = 35, h_to_unit :: String = "ft", V_2_coef = 1.2, μ = NaN)
    T = promote_type(eltype(aircraft.N_e), eltype(aircraft.engine.λ_bpr), eltype(aircraft.engine.prop_disc_load), eltype(aircraft.engine.k_c), eltype(aircraft.CL_max), eltype(TODA), eltype(α), eltype(β), eltype(h), eltype(γ_2_min))
    @assert h >= 0 && h_to >= 0 && V_2_coef >= 0  && TODA >= 0 "h, V and takeoff distance inputs must be non-negative"
    @assert all(>=(0), any(isnan, α)) && all(<=(1), any(isnan, α)) && β >= 0 && β <= 1 "α and β must be between 0 and 1"
    
    # If no friction value was specified
    if isnan(μ)
        # Get the runway data
        runway = ParseData.get_data(data = "Runway")

        # Tries to find the surface data inside the dataset, if cannot find, 
        # then throws an error!
        μ = try
            runway[only(findall(==(surface), runway.Surface)), "Friction"]
        catch e
            throw("$surface data cannot be found!")
        end
    else
        # If a value for friction was added, check its validity
        @assert μ >= 0 "μ input must be non-negative"
    end

    return TakeOffPerformance{T}(aircraft, TODA, TODA_unit, α, β, h, h_unit, γ_2_min, γ_2, h_to, h_to_unit, V_2_coef, μ)
end

# Takeoff performance calculation: 
# WS_min = Minimum Wing Loading considered
# WS_max = Maximum Wing Loading considered
# WS_N = Number of Wing Loading points
function bfl_perf(
    to_data :: TakeOffPerformance,
    WS_min = 1,
    WS_max = 10000,
    WS_N = 1000,
    )

    WS = LinRange(WS_min*u"kg/m/s^2", WS_max*u"kg/m/s^2", WS_N)
    WS_imp = uconvert.(u"lbf/ft^2", WS)

    ### Convert units
    # Convert to unit form
    h_unit = AeroUnits.convert_to_unit(to_data.h_unit)
    h = uconvert(u"ft", to_data.h * h_unit)
    h_to_unit = AeroUnits.convert_to_unit(to_data.h_to_unit)
    h_to = uconvert(u"ft", to_data.h_to * h_to_unit)
    TODA_unit = AeroUnits.convert_to_unit(to_data.TODA_unit)
    BFL = uconvert(u"ft", to_data.TODA * TODA_unit)

    ### Determine density
    (ρ,_,_,_) = ISAdata(h)
    (ρ_0,_,_,_) = ISAdata(0u"ft")

    σ = ρ / ρ_0
    ρ = uconvert(u"slug/ft^3", ρ)

    ### Determine whether the input is power or thrust 
    if to_data.aircraft.engine.type == "Jet"
        k_e = 0.75 * ((5 + to_data.aircraft.engine.λ_bpr) / (4 + to_data.aircraft.engine.λ_bpr))
    else
        prop_disc_load = to_data.aircraft.engine.prop_disc_load * u"hp/ft^2"
        k_e = 5.75 * to_data.aircraft.engine.k_c * cbrt(σ / prop_disc_load)
    end

    g = uconvert(u"ft/s^2", 1*u"ge") # Gravitational acceleration constant
    α = to_data.α
    β = to_data.β
    CL_max_to = to_data.aircraft.CL_max + to_data.aircraft.ΔCL_max_to
    CL_2 = CL_max_to / (to_data.V_2_coef ^ 2)
    Δγ2 = to_data.γ_2 - to_data.γ_2_min
    μ_prime = 0.01 * to_data.aircraft.CL_max + to_data.μ

    TW = (α ./ (β * k_e)) .* ((1 ./ ((1.159 * (BFL - 655u"ft" / sqrt(σ)) * (1 + 2.3 * Δγ2)) ./ (((WS_imp .* α) ./ (ρ * g * CL_2) .+ h_to)) .- 2.7)) .+ μ_prime)

    if to_data.aircraft.engine.type == "Propeller"
        TW = upreferred.(TW)
        TW = ustrip.(TW) * u"m/s"
    end
    return (WS, TW)
end

##### Landing Performance Definition
abstract type AbstractLandingPerformance end

struct LandingPerformance{T} <: AbstractLandingPerformance
    aircraft    :: Aircraft{T}
    LDA         :: T # LDA = Landing Distance Available (default LDA_unit = ft)
    LDA_unit    :: String
    α           :: T # α = current weight / total weight ratio
    h           :: T # h = height of aircraft (default h_unit = ft)
    h_unit      :: String
    S_a         :: T # Obstacle Clearance distance
    S_a_unit    :: String
    V_land_coef :: T # V_land_coef = Landing speed as a function of stall speed, usually 1.15 for FAR and 1.1 for military
    K_r         :: T # Deceleration factor
    K_l         :: T # Multiplier for landing distance (3/5 for FAR25)
end

function LandingPerformance(; aircraft :: Aircraft, LDA, LDA_unit :: String = "ft", α, h = 0, h_unit :: String = "ft", S_a = 305, S_a_unit :: String = "m", V_land_coef = 1.15, K_r = 1, K_l = 3/5)
    T = promote_type(eltype(aircraft.CL_max), eltype(LDA), eltype(α), eltype(h), eltype(S_a), eltype(V_land_coef), eltype(K_r), eltype(K_l))
    @assert h >= 0 && S_a >= 0 && LDA >= 0 "h, S_a and landing distance inputs must be non-negative"
    @assert α >= 0 && α <= 1 "α must be between 0 and 1"

    return LandingPerformance{T}(aircraft, LDA, LDA_unit, α, h, h_unit, S_a, S_a_unit, V_land_coef, K_r, K_l)
end

# Landing performance calculation: 
function landing_perf(landing_data :: LandingPerformance)
    ### Convert units
    # Convert to unit form
    h_unit = AeroUnits.convert_to_unit(landing_data.h_unit)
    h = uconvert(u"m", landing_data.h * h_unit)
    LDA_unit = AeroUnits.convert_to_unit(landing_data.LDA_unit)
    LDA = uconvert(u"m", landing_data.LDA * LDA_unit)
    S_a_unit = AeroUnits.convert_to_unit(landing_data.S_a_unit)
    S_a = uconvert(u"m", landing_data.S_a * S_a_unit)

    ### Determine density
    (ρ,_,_,_) = ISAdata(h)
    (ρ_0,_,_,_) = ISAdata(0u"m")

    σ = ρ / ρ_0

    g = uconvert(u"m/s^2", 1*u"ge") # Gravitational acceleration constant
    α = landing_data.α
    K_l = landing_data.K_l
    K_r = landing_data.K_r
    CL_max_land = landing_data.aircraft.CL_max + landing_data.aircraft.ΔCL_max_ld

    WS = (((LDA * K_l - S_a)*(σ * CL_max_land)) ./ ((5/g) .* K_r .* α))
    WS = ustrip(WS) * u"kg / m / s^2"

    return WS
end

function update_constraint_region(;WS_new, boundary = [], TW_new = 0, WS_max = 999999*u"kg/m/s^2")
    if boundary == []
        boundary = DataFrame(WS = WS_new, TW_bound = TW_new)
    elseif TW_new == 0
        WS_max = min(WS_max, WS_new)
    else
        df_new = DataFrame(WS = WS_new, TW = TW_new)
        boundary = outerjoin(boundary, df_new, on = :WS)
        boundary.TW_bound = max.(boundary.TW_bound, boundary.TW)
        boundary = boundary[:, ["WS", "TW_bound"]]
    end

    boundary = boundary[boundary.WS .<= WS_max, :]

    return (boundary, WS_max)
end

function create_constraint_diagram(;p = nothing, aircraft :: Aircraft, mission :: MissionProfile, V_cruise, W_frac, α_land, WS_start, WS_end, WS_num, field_condition, TW, total_cost, show_cost :: Bool = false)
    WS = LinRange(WS_start, WS_end, WS_num)
    
    # For plotting
    if isnothing(p)
        draw_constraint = false
    else
        draw_constraint = true
    end

    # Initialise key variables to save information
	global stage_num = 1;
	global α_running = ones(length(WS_num),1);

	global α_save = Array{Float64}(undef, WS_num, length(mission.stage));
	global WS_max_save = WS_end*u"kg/m/s^2";
	global boundary_save = [];
	global cruise_idx2 = 1;

    if show_cost == true && draw_constraint == true
	    Plots.contourf!(p, WS.*u"kg/m/s^2", TW, total_cost', color = :turbo, levels = 20)
    end

	# For loop for each stage specified
	for stage in mission.stage
		if stage == "Takeoff"
			takeoff_var = ConstraintDiagram.TakeOffPerformance(
			    aircraft   = aircraft,
			    TODA       = mission.distance[stage_num],
			    TODA_unit  = "m",
			    α          = α_running,
			    β          = 1,
			    h          = mission.h[stage_num],
				h_unit     = "m",
			    surface    = field_condition[1]
			)

			(WS_plot, TW_plot) = ConstraintDiagram.bfl_perf(takeoff_var, WS_start, WS_end, WS_num);
		elseif stage == "Landing"
			if field_condition[2] == false
				K_r = 1
			else
				K_r = 0.66 # If reverser is allowed
			end
			
			landing_var = ConstraintDiagram.LandingPerformance(
			    aircraft   = aircraft,
			    LDA        = mission.distance[stage_num],
			    LDA_unit   = "m",
			    α          = α_land,
			    h          = mission.h[stage_num],
				h_unit     = "m",
			    K_r        = K_r
			)

			WS_plot = ConstraintDiagram.landing_perf(landing_var);
			#global WS_land = WS_plot
		elseif stage == "Cruise"
			point_var = ConstraintDiagram.PointPerformance(
			    aircraft   = aircraft,
			    α          = α_running,
			    β          = 1,
			    h          = mission.h[stage_num],
				h_unit     = "m",
			    V∞         = V_cruise[:,cruise_idx2],
				V∞_unit    = "m/s",
				G          = mission.G[stage_num],
				N_a        = mission.N_a[stage_num]
			)

			(WS_plot, TW_plot) = ConstraintDiagram.point_perf(point_var, WS_start, WS_end, WS_num);
		else
			point_var = ConstraintDiagram.PointPerformance(
			    aircraft   = aircraft,
			    α          = α_running,
			    β          = 1,
			    h          = mission.h[stage_num],
				h_unit     = "m",
			    V∞         = [mission.V∞[stage_num]],
				V∞_unit    = "m/s",
				G          = mission.G[stage_num],
				N_a        = mission.N_a[stage_num]
			)

			(WS_plot, TW_plot) = ConstraintDiagram.point_perf(point_var, WS_start, WS_end, WS_num);
		end

		# Landing is a specific case where it is a vertical plot
        if stage == "Landing"
            if draw_constraint == true
                vline!(p, [WS_plot], label = mission.name[stage_num])
            end
            (boundary, WS_max) = ConstraintDiagram.update_constraint_region(boundary = boundary_save, WS_new = WS_plot, WS_max = WS_max_save)
        else
            if draw_constraint == true
                Plots.plot!(p, WS_plot, TW_plot, label = mission.name[stage_num])
            end
            (boundary, WS_max) = ConstraintDiagram.update_constraint_region(boundary = boundary_save, WS_new = WS_plot, TW_new = TW_plot[:,1], WS_max = WS_max_save)
        end

		# Update global variables
		if stage == "Cruise"
			global α_running = α_running .* W_frac[:,cruise_idx2]
			global cruise_idx2 = cruise_idx2 + 1
		else
			global α_running = α_running * mission.w_frac[stage_num]
		end
		
		α_save[:, stage_num] .= α_running
		
		global stage_num = stage_num + 1
		global boundary_save = boundary
		global WS_max_save = WS_max
	end

    if length(WS) == 1
        cost_interp = linear_interpolation(TW, total_cost[1,:], extrapolation_bc = +Inf)
        boundary_cost = cost_interp(ustrip.(boundary_save.TW_bound))
    else
        cost_interp = linear_interpolation((WS, TW), total_cost, extrapolation_bc = +Inf)
        boundary_cost = cost_interp(ustrip.(boundary_save.WS), ustrip.(boundary_save.TW_bound))
    end

    boundary_cost[isnan.(boundary_cost)] .= +Inf

    if findmin(boundary_cost)[1] == +Inf
        min_cost = NaN
        WS_des = NaN
        TW_des = NaN
    else
        idx_min = findmin(boundary_cost)
        min_cost = idx_min[1]
        if min_cost < 0 # TODO: Understand why this error occur
            min_cost = NaN
            WS_des = NaN
            TW_des = NaN
        else
            WS_des = ustrip.(boundary_save.WS[idx_min[2][1]])
            TW_des = ustrip.(boundary_save.TW_bound[idx_min[2][1]])
        end
    end

    if draw_constraint == true
        # Specify title
        title!(p, "Constraint Diagram")
        xlabel!(p, join(["Wing Loading (W/S) ", u"kg / m / s ^ 2"]))

        # Change the label depending on engine type
        if aircraft.engine.type == "Jet"
            ylabel!(p, "Thrust to Weight Ratio (T/W)")
            ylims!(p, 0, 1)
        else
            ylabel!(p, join(["Power to Weight Ratio (P/W) ", u"m / s"]))
            ylims!(p, 0, 500)
        end

        Plots.plot!(p, ustrip.(boundary_save.WS), ustrip.(boundary_save.TW_bound),
            linestyle=:dash, linewidth=5, label = "Boundary")

        Plots.scatter!(p, ustrip.([WS_des]), ustrip.([TW_des]), markersize = 5, label = "Design Point")
    else
        return (WS_des, TW_des, min_cost)
    end
end

end
