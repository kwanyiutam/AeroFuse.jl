# Main module for cost model
module CostModel

# Initialise packages used
using ISAData
using Unitful
using CSV
using DataFrames

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the ParseData module
include("parsedata.jl")

# Include the InputValidate module
include("inputvalidate.jl")

"""
    `cost_init` - A function which initiates the dataframe with the cost model

    Returns the cost model in a dataframe form
"""
function cost_init()
    # Get cost and inflation data
    df_cost = ParseData.get_data(data = "Cost")
    df_inflation = ParseData.get_data(data = "Inflation")

    # Divide the maximum CPI value (from most recent year) with the CPI of that years
    max_CPI = df_inflation[argmax(skipmissing(df_inflation[!, :Year])), :CPI]
    df_inflation[:, :CPI] = max_CPI ./ df_inflation[:, :CPI]

    # Rename the CPI column to factor
    rename!(df_inflation, :CPI => :Factor)

    # Convert the column to String so that it will be able to join
    df_inflation[!, :Year] = Vector{String}(string.(df_inflation[!, :Year]))

    # Join the two dataframes together
    leftjoin!(df_cost, df_inflation, on = :Year)

    return df_cost
end

function cost_model(;df_cost::DataFrame, model::String = "Full")
    # Filter the cost functions by model
    row_idx = findall(value -> occursin(Regex("(?i)"*model),value),df_cost[:,"Model"])

    # Throw error if the model cannto be found
    if isempty(row_idx)
        throw(ArgumentError("The model chosen ($model) cannot be found on the cost model!"))
    end

    # Get the row with Units
    unit_row_idx = findfirst(==("Units"),df_cost[:,1])

    # Concatenated the two rows
    row_idx = vcat(unit_row_idx,row_idx)

    # Get a cleaned copy, without other costs
    df_clean = deepcopy(df_cost[row_idx,:])

    # Drop any columns with all missing data, except for the units row
    cols_to_drop = names(df_clean)[all.(ismissing, eachcol(df_clean[setdiff(1:end, unit_row_idx), :]))]
    select!(df_clean, Not(cols_to_drop))

    return df_clean
end

end


