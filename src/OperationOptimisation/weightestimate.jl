module WeightEst

# Initialise packages used
using DataFrames
using Unitful

# Include the InputValidate module
include("inputvalidate.jl")

"""
    `weight_calculation` - A function which calculates a weight estimate

    return the weight estimation in dataframe
"""
function weight_calculation(;df_aircraft::DataFrame,N_aircraft::Int,aircraft_idx::Int,df_mission::DataFrame,N_stages::Int,df_payload::DataFrame,N_payload::Int,payload_idx::Int)



    return df_aircraft
end

end