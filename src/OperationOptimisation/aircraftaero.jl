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

    return df_aircraft
end

end