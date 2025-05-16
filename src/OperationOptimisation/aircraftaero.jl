module AircraftAero

# Initialise packages used
using DataFrames
using Unitful
using AeroFuse

# Include the InputValidate module
include("inputvalidate.jl")

"""
    `run_aero_analysis` - A function which runs the aerodynamics analysis for OperationOptimisation

    return the original dataframe
"""
function run_aero_analysis(;df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int)
    #taper_ratio = df_aircraft[]
    ## Initial guess
    n_vars = 64 # Number of spanwise stations
    c = 0.125 # Fixed chord
    c_w = LinRange(c, c, n_vars) # Constant distribution
    CL_tgt = 1.6  # Target lift coefficient
    nc = length(c_w)

    return df_aircraft
end

end