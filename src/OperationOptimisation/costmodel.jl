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

function cost_variables(;variables::Vector{String}, units, df_aircraft::DataFrame, aircraft_idx::Int, df_mission::DataFrame, N_stages::Int, df_payload::DataFrame, payload_idx::Int)
    if length(variables) != length(units)
        throw(ArgumentError("Variables and Units vector do not match! Please check..."))
    end

    values = []

    for var_idx in eachindex(variables)
        var = variables[var_idx]

        if var == "Fuel Weight"
            value = df_aircraft[findfirst(==("Fuel Weight"),df_aircraft[:,1]),aircraft_idx]
        elseif var == "Cruise Time"
            value = sum(df_mission[findfirst(==("Duration"),df_mission[:,1]),ncol(df_mission)-N_stages+1:ncol(df_mission)])
        elseif var == "Number of Cabin Crew"
            value = df_payload[findfirst(==("Cabin Crew"),df_payload[:,1]),payload_idx]
        elseif var == "MTOW"
            value = df_aircraft[findfirst(==("MTOW"),df_aircraft[:,1]),aircraft_idx]
        else
            throw(ArgumentError("Invalid Cost Variable $var ! This variable has not been defined in the cost model, please update the cost_variables function!"))
        end

        unit = units[var_idx]

        if !ismissing(unit)
            value = uconvert(unit, value)
        end

        push!(values, value)
    end

    return values
end

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

    constant_col = findfirst(==("Constant"),DataFrames.names(df_cost))

    # Convert all the parameters into any vector so that they can be parsed
    for col in constant_col+1:ncol(df_cost)-1
        df_cost[!, col] = convert(Vector{Any}, df_cost[!, col])
    end

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

    unit_row_again = findfirst(==("Units"),df_clean[:,1])
    constant_col = findfirst(==("Constant"),DataFrames.names(df_clean))

    for col in constant_col+1:ncol(df_clean)
        unit = df_clean[unit_row_again,col]
        if !ismissing(unit)
            df_clean[unit_row_again,col] = AeroUnits.convert_to_unit(unit)
        end
    end

    return df_clean
end

function cost_calc(;df_cost_model::DataFrame, df_aircraft::DataFrame, aircraft_idx::Int, df_mission::DataFrame, N_stages::Int, df_payload::DataFrame, payload_idx::Int)
    # Get the row with Units
    unit_row_idx = findfirst(==("Units"),df_cost_model[:,1])

    # Rows which are cost calculations
    cost_rows = setdiff(1:nrow(df_cost_model), unit_row_idx)

    # Define a new dataframe to save information
    df_results = DataFrame(Cost_Component = df_cost_model[cost_rows,1], Category = df_cost_model[cost_rows,"Category"])
    df_results[!, "Cost (USD)"] = fill(0.0, nrow(df_results))

    # Get the constant column and variables
    constant_col_idx = findfirst(==("Constant"), DataFrames.names(df_cost_model))
    variable_col_idx = constant_col_idx+1:(ncol(df_cost_model)-1)
    variables = DataFrames.names(df_cost_model)[variable_col_idx] # Factor column is at the back
    units = collect(df_cost_model[unit_row_idx,variable_col_idx])

    values = cost_variables(variables = variables, units = units, df_aircraft=df_aircraft, aircraft_idx=aircraft_idx, df_mission=df_mission, N_stages=N_stages, df_payload=df_payload, payload_idx=payload_idx)
    
    # Replace missing values with zero
    values = ustrip.(coalesce.(values, 0.0))

    for idx in cost_rows
        component = df_cost_model[idx, 1]
        logic = df_cost_model[idx, "Logic"]
        
        # Replace missing values with zero
        coefficients = collect(df_cost_model[idx, variable_col_idx])
        coefficients = map(x -> ismissing(x) ? 0.0 : parse(Float64, string(x)), coefficients)

        if logic == "SUM"
            # Dot product then summation to calculate the cost, plus constant
            cost = sum(values' * coefficients) + df_cost_model[idx, constant_col_idx]
        elseif logic == "PROD"
            # Constant multiplier
            cost = df_cost_model[idx, constant_col_idx]

            for idx in eachindex(coefficients)
                cost *= values[idx]^(coefficients[idx]) 
            end
        elseif logic == "FCREW"
            MTOW = values[findfirst(==("MTOW"), variables)]
            time = values[findfirst(==("Cruise Time"), variables)]
            cost = (0.000326*MTOW+653)*time
        else
            throw(ArgumentError("Logic cannot be parsed, please double check!"))
        end

        df_row_idx = findfirst(==(component), df_results[:,1])
        df_results[df_row_idx,"Cost (USD)"] = cost * df_cost_model[idx, "Factor"]
    end

    total_cost = sum(collect(df_results[:,"Cost (USD)"]))

    return (total_cost, df_results)
end

end