module InputValidate

# Initialise packages used
using CSV
using DataFrames
using AeroFuse
using Unitful

# Initialise the aero units
include("aerounits.jl")
Unitful.register(AeroUnits)

# Include the ParseData module
include("parsedata.jl")

function get_value(df::DataFrame,label,col=nothing,all_col::Bool=false)
    if all_col == true
        value = df[findfirst(==(label),df[:,1]),:]
    else
        value = df[findfirst(==(label),df[:,1]),col]
    end

    return value
end

"""
    `rowcol_names_unique` - A function which checks whether the rows or columns names were unique
    
    If not, it will throw an exception to ask for the user to check the inputs
"""
function rowcol_names_unique(df::DataFrame)
    if(length(names(df)) != length(unique(names(df))))
        throw(ErrorException("Column names have been duplicated! Please check your inputs!"))
    end

    if((length(df[:,1])) != length(unique(df[:,1])))
        throw(ErrorException("Row names have been duplicated! Please check your inputs!"))
    end
end

"""
    `change_column_type` - A function which changes all the columns into a string+missing column
    
    This ensure that the columns allow missing and string values to coexist, because sometimes this is allowed
"""
function change_column_type(df)
    # For each column in the dataframe
    for col in DataFrames.names(df)
        # If the type is InlineString, this is too unstable
        # If the type is Missing, then it needs to accept string!
        if (eltype(df[!, col]) <: InlineString) || (eltype(df[!, col]) <: Missing)
            df[!, col] = convert(Vector{Union{Missing, String}}, df[!, col])
        end
    end

    return df
end

"""
    `rowcolumn_check` - A function which checks whether the rows or columns contains relevant information

    Takes `df` as the dataframe to check, and `df_assumptions` as the reference dataframe.

    `array_check` is an array consisting of the variables which should exist as row/columns names.

    `mode` can be used to specify further action if the data cannot be found.

    `append_rowcol` is a string to specify which row or column of data to append (will search for that word)

    return the validation checks
"""
function rowcolumn_check(; df::DataFrame, df_assumptions::DataFrame, array_check::Vector{String}, axis::String = "Row", mode::String = "Mandatory", append_rowcol::String = "Append", verbose::Bool = false)
    # Check that the axis inputs are correct
    if !(axis in ["Row", "Column"])
        throw(ArgumentError("Axis must be 'Row' or 'Column'"))
    end

    # If the input is row, check = first column, else, it is the column names
    check = axis == "Row" ? df[:, 1] : names(df)

    # For each check in the array
    for i in array_check
        # Find the first instance of the value occuring
        idx = findfirst(==(i), check)

        # If nothing can be found
        if isnothing(idx)
            # If the data is mandatory, throw an exception
            if mode == "Mandatory"
                throw(ErrorException("$i cannot be found in the $axis !"))
            # If the data can be appended with defaults
            elseif mode == "Append"
                # And if it is for the rows
                if axis == "Row"
                    # Find the first instance of the array_check value occuring in df_assumptions
                    row_idx = findfirst(==(i), df_assumptions[:, 1])

                    if verbose == true
                        @warn "$i not found in input, adding default values from sample."
                    end

                    # Change the dataframe to ensure that it can take string and missing values
                    df = change_column_type(df)

                    # Take the first row (as a template)
                    row = deepcopy(df[1, :])

                    # Append corresponding values in each colum into the row
                    for col in DataFrames.names(row)
                        # Find the matched column name (e.g., "Units")
                        idx = findfirst(==(col), DataFrames.names(df_assumptions))

                        # If the column name cannot be found, this is assumed to be
                        # User specified configurations
                        if idx === nothing
                            # Append the information from the append_rowcol column instead
                            idx = findfirst(==(append_rowcol), DataFrames.names(df_assumptions))

                            # If still cannot find, something is wrong with append_rowcol
                            if idx === nothing
                                throw(ArgumentError("Invalid inputs for append_rowcol, cannot find $append_rowcol !"))
                            end
                        end
    
                        # Add the data into the row
                        row[col] = df_assumptions[row_idx, idx]
                    end

                    # Append the data to the back
                    push!(df, row, promote = true)
                else
                    col_idx = findfirst(==(i), DataFrames.names(df_assumptions))

                    if verbose == true
                        @warn "$i not found in input, adding default values from sample."
                    end

                    # Copy the column for a foundation
                    col = deepcopy(df_assumptions[:, col_idx])
                    names(col) = i

                    # Convert the column type to be the same as the one in df_assumptions
                    #col = convert(Vector{eltype(df_assumptions[!, i])}, col)

                    # Append corresponding values in each row into the column
                    for row in 1:nrow(df)
                        # Find the matched row name (e.g., "Stage")
                        idx = findfirst(==(df[row, 1]), df_assumptions[:, 1])

                        # If the Row name cannot be found, this is assumed to be
                        # Extra rows specified by the user
                        if idx === nothing
                            # Append the information from the append_rowcol column instead
                            idx = findfirst(==(append_rowcol), df_assumptions[:, 1])

                            # If still cannot find, something is wrong with append_rowcol
                            if idx === nothing
                                name = df[row, 1]
                                throw(ArgumentError("Invalid inptus for append_rowcol, cannot find $name !"))
                            end
                        end
    
                        # Add the data into the column
                        col[row] = df_assumptions[idx, col_idx]
                    end

                    df[!, i] = col
                end
            end
        end
    end

    return df
end

"""
    `replace_reference` - A function which replaces references with values

    Takes `df` as the dataframe to check, and `key_cols` as the columns NOT to be included.

    For example, if two aircraft, "Baseline" and "Alternate" were defined, and the user
    defined one of the values as "Baseline", then replace that with the "Baseline"'s value in
    that row.s

    return the replaced dataframe
"""
function replace_reference(;df::DataFrame,key_cols::Vector{String},data_type::String,skip_row::Vector{Int}=Int[],verbose::Bool = false)
    # Get all configuration names
    data_idx = findall(!in(key_cols), DataFrames.names(df))
    data_configs = DataFrames.names(df)[data_idx]

    # For each row in the dataframe and each data
    for row in 1:nrow(df)
        # If row needs to be skipped, then skip
        if row in skip_row
            continue
        end
    
        #For each column
        for idx in data_idx
            # Get the value of the cell
            value = df[row,idx]

            # If the data is missing, skip
            if ismissing(value)
                continue
            end

            # See if the input matches with any of the aircradataft names
            match_data = findfirst(==(value), data_configs)

            # If matched
            if !isnothing(match_data)
                # Get the matched data index in the df
                match_data = data_idx[match_data]

                # If the match_data index is the same as the data
                if match_data == idx
                    throw(ErrorException("$data_type '$value' input in Row $row is referring to themselves, which cannot happen (no self-referencing), please check your inputs!"))
                else
                    # Otherwise, append data from the matched config
                    name_current = DataFrames.names(df)[idx]

                    if verbose == true
                        @info "Appending $data_type '$value' data for row $row into $data_type '$name_current'"
                    end

                    # Add the data from the matched data
                    df[row,idx] = df[row,match_data]
                end
            end
        end
    end

    return df
end

"""
    `operator_compare` - A function which checks mathematical operators

    Takes `df` as the dataframe to check, `rule` as the given rule, `operator` as the known operator,
    `row` as the row number of the data and `N_config` as the number of configurations specified (i.e, extra columns)

    no return, it is a checking function (void)
"""
function operator_compare(;df::DataFrame,rule::String,operator::String,row::Int,N_config::Int)
    # Checks to see whether the inputs are sensisble
    if N_config < 1
        throw(ArgumentError("Number of configurations specified must be bigger than 1!"))
    end

    value = parse(Float64,rule[length(operator)+1:end])

    for col in (ncol(df)-N_config+1):ncol(df)
        data = df[row,col]
        if ismissing(data)
            continue
        end

        data = parse(Float64,df[row,col])

        assumption = df[row,"Assumptions"]
        config = DataFrames.names(df)[col]

        if operator == ">="
            @assert data >= value "Value at Row $assumption Column $config must be bigger than or equal to $value"
        elseif operator == "<="
            @assert data <= value "Value at Row $assumption Column $config must be less than or equal to $value"
        elseif operator == ">"
            @assert data > value "Value at Row $assumption Column $config must be bigger than $value"
        elseif operator == "<"
            @assert data < value "Value at Row $assumption Column $config must be less than $value"
        elseif operator == "=="
            @assert data == value "Value at Row $assumption Column $config must be equal to $value"
        elseif operator == "!="
            @assert data != value "Value at Row $assumption Column $config must not be equal to $value"
        else
            throw(ArgumentError("Cannot parse the operator at Row $assumption Column $config"))
        end
    end
end

"""
    `design_check` - A function which defines the design parameters

    Takes `df` as the dataframe to check for design parameters, and `design_param` to save them
    
    `location` specifies where the design parameters originated from (to trace the origin)

    `N_config` defines where the beginning of the columns to check

    Returns the two dataframes
"""
function design_check(;df::DataFrame,design_param::DataFrame,location::String,N_config::Int)
    # Checks to see whether the inputs are sensisble
    if N_config < 1
        throw(ArgumentError("Number of configurations specified must be bigger than 1!"))
    end

    # Specify types of design parameters, optimisation or design
    types_regex = ["OPT;","DES;","PERT;","MIN_COST;","MIN_FUEL"]
    types = ["Optimise", "Design", "Perturbations", "Minimum Cost Velocity", "Minimum Fuel Velocity"]

    # For each column
    for col in (ncol(df)-N_config+1):ncol(df)
        # For each design type
        for idx in eachindex(types_regex)
            # Check whether the matching regex has been found
            check = findall(value -> !ismissing(value) && occursin(Regex(types_regex[idx]),value),df[:,col])
            
            # If some data was found
            if !isempty(check)
                # For each check
                for i in check
                    # Get the data and split by ;
                    data = df[i,col]
                    data = split(data,";")

                    # Get the upper and lower boundary
                    lb = parse(Float64,data[3])
                    ub = parse(Float64,data[4])

                    # Save the design, following the order: Design Parameter => Saved Location
                    # => Saved Row => Saved Column => Type => Continuous => Lower Bound => Upper Bound
                    save_design = [data[2],location,i,col,types[idx],df[i,"Variable Type"]=="Float64",lb,ub]
                    push!(design_param,save_design)

                    # Replace the assumption table with lower bound, just as a placeholder for now
                    df[i,col] = data[3]
                end
            end
        end
    end

    return (df, design_param)
end

"""
    `mandatory_check` - A function which checks whether all mandatory variables were specified

    Takes `df` as the dataframe to check, and `key_col` as the column name with the mandatory column

    no return, it is a checking function (void)
"""
function mandatory_check(;df::DataFrame,key_col::String = "Logic",N_config::Int)
    # Checks to see whether the inputs are sensisble
    if N_config < 1
        throw(ArgumentError("Number of configurations specified must be bigger than 1!"))
    end

    # Convert to boolean column
    df[!, key_col] = convert(Vector{Bool}, df[!, key_col])

    # Find all true indicies and get those rows from the specified columns
    idx_true = findall(==(true),df[!, key_col])    
    df_check = df[idx_true, (ncol(df)-N_config+1):ncol(df)]

    # Check whether any values are missing
    check = any(ismissing, Iterators.flatten(eachcol(df_check)))

    # If missing, throws error
    if check == true
        throw(ArgumentError("A mandatory variable is missing! Please double check your inputs."))
    end
end
    

"""
    `parse_validation` - A function which conducts logic checks for the inputs

    Takes `df` as the dataframe to check, and `key_col` as the column name with the logic rules

    Return adapted and validated date
"""
function parse_validation(;df::DataFrame,key_col::String = "Logic",N_config::Int)
    # Checks to see whether the inputs are sensisble
    if N_config < 1
        throw(ArgumentError("Number of configurations specified must be bigger than 1!"))
    end

    # For each row in the dataframe
    for row in 1:nrow(df)
        # Check the rules in that row
        rules = df[row,key_col]

        # If the data is missing, then skip (no rules)
        if ismissing(rules)
            continue
        end

        # Else, split the rules up by ";" delimiter
        rules = split(rules, ";")
        for rule in rules
            # Clean the rules to remove any spaces
            rule = strip(rule,' ')

            # If the rule is airfoildata
            if lowercase(rule) == "airfoildata"
                # Get the list of airfoils in the database
                airfoils_list = readdir("./airfoil_database")

                # For each defined column
                for col in (ncol(df)-N_config+1):ncol(df)
                    # Get the data
                    data = df[row,col]

                    # If the data is a NACA airfoil 4 digit OR exist in our database, continue
                    if occursin(r"(?i)^NACA\d{4}$", data)
                        continue
                    elseif any(s -> occursin(Regex("(?i)^"*data*".dat"), s), airfoils_list)
                        continue
                    end

                    assumption = df[row,"Assumptions"]
                    config = DataFrames.names(df)[col]
                    
                    # Throw an error if the airfoil cannot be found
                    throw(ErrorException("Airfoil $data defined at Row $assumption Column $config cannot be found in our database!"))
                end
            # If the data started with Data (Get data from other sources)
            elseif occursin(r"(?i)^Data:", rule)
                # Get the reference value
                reference = split(rule,":")[2]
                df_ref = ParseData.get_data(data = String(reference))

                # First three must be the assumptions, units and variable types
                ref_data = df_ref[:,["Assumptions","Units","Variable Type"]]

                # For each defined data
                for col in (ncol(df)-N_config+1):ncol(df)
                    # Get the reference value and match it in the data
                    ref_value = df[row,col]
                    idx = findfirst(==(ref_value),DataFrames.names(df_ref))

                    # Throw error if the reference value cannot be found
                    if idx === nothing
                        assumption = df[row,"Assumptions"]
                        config = DataFrames.names(df)[col]
                        
                        throw(ArgumentError("Reference value $ref_value in Row $assumption Column $config cannot be found in Dataset $reference"))
                    end
                    
                    # Append the column into the daataset
                    ref_data[!, DataFrames.names(df)[col]] = df_ref[:,idx]
                end

                # Fill the missing data for category, logic and mandatory
                ref_data[!, "Category"] .= reference
                ref_data[!, "Logic"] .= missing
                ref_data[!, "Mandatory"] .= true # Placeholder at this point, not required

                # Append the dataset into the assumptions
                append!(df, ref_data, promote = true)
            # If the data started with Match (Inputs should match the options)
            elseif occursin(r"(?i)^Match:", rule)
                # Get the regex inside the double quotes
                regex_match = match(r"\"([^\"]+)\"", rule)

                # If the regex_match is not found, then the rule is in the wrong format
                if regex_match === nothing
                    throw(ArgumentError("Invalid format: '$rule'"))
                end

                # Get the first item of the regex match
                regex_match = regex_match[1]

                # For each defined column
                for col in (ncol(df)-N_config+1):ncol(df)
                    data = df[row,col]

                    # If the data is missing, then it is not acceptable
                    if ismissing(data)
                        assumption = df[row,"Assumptions"]
                        config = DataFrames.names(df)[col]
                        
                        throw(ArgumentError("Data for Row $assumption Column $config is missing, please check!"))
                    elseif !occursin(Regex(regex_match), data)
                        assumption = df[row,"Assumptions"]
                        config = DataFrames.names(df)[col]
                        
                        throw(ArgumentError("Value '$data' at Row $assumption Column $config does not match the pattern $regex_match"))
                    end
                end
            # If the data started with Req (Data required)
            elseif occursin(r"(?i)^Req:", rule)
                # Get the regex inside the brackets, separated by a comma
                requirement = match(r"Req:\(\"([^\"]+)\",\"([^\"]+)\"\)", rule)

                # If the requirements is not found, then the rule is in the wrong format
                if requirement === nothing
                    throw(ArgumentError("Invalid format: '$rule'"))
                end

                # Captrues the label and regex from the bracket
                label, regex = requirement.captures

                idx = findfirst(==(label),df[:,1])

                # For each defined column
                for col in (ncol(df)-N_config+1):ncol(df)
                    data = df[idx,col]

                    # If the regex matched (i.e. an input is required) but the input is also missing, error
                    if occursin(Regex(regex),data) && ismissing(df[row,col])
                        assumption = df[row,"Assumptions"]
                        config = DataFrames.names(df)[col]
                        
                        throw(ArgumentError("Value at Row $assumption Column $config is missing, but is required since $label '$data' satisfies $regex"))
                    end
                end
            # If this is a known mathematical operator
            elseif occursin(r"^(>=|<=|>|<|==|!=)-?\d+", rule)
                # Get the operator
                get_operator = match(r"^(>=|<=|>|<|==|!=)", rule)[1]
                # Check whether the values satisfies the operators
                operator_compare(df = df,rule = String(rule),operator = String(get_operator), row = row, N_config = N_config)
            else
                throw(ErrorException("Logic $rule cannot be understood on Row $row, please check the logic or submit a request for additional rules"))
            end
        end
    end

    return df
end

"""
    `mission_check` - A function which checks mission inputs.

    Checks order of the missions specified AND the mission-specific inputs

    no return, it is a checking function (void)
"""
function mission_check(; df :: DataFrame, N_config::Int)
    # Checks to see whether the inputs are sensisble
    if N_config < 1
        throw(ArgumentError("Number of configurations specified must be bigger than 1!"))
    end

    for col in (ncol(df)-N_config+1):ncol(df)
        stage_row = findfirst(==("Stage"),df[:,"Assumptions"])
        gradient_row = findfirst(==("Gradient"),df[:,"Assumptions"])
        stage = lowercase(df[stage_row,col])
        order = col - (ncol(df)-N_config+1)

        if (order == 0) && (stage != "takeoff")
            throw(ArgumentError("Invalid First Stage: Stage was specified as $stage but it must start with takeoff"))
        end

        if (order == N_config-1) && (stage != "landing")
            throw(ArgumentError("Invalid Last Stage: Stage was specified as $stage but it must end with landing"))
        end

        # Get the previous stage
        prev_stage = lowercase(df[stage_row,col-1])

        # Checks different conditions to see if the order and values make sense
        if stage == "takeoff"
            @assert (prev_stage == "landing") || (order == 0) "If takeoff is not the first stage, it must be preceded by landing!"
        elseif stage == "climb"
            @assert (prev_stage == "takeoff" || prev_stage == "descend") "Climb must be preceded by takeoff or descend stage!"
            @assert (parse(Float64,df[gradient_row,col]) > 0) "Climb gradient must be bigger than zero!"
        elseif stage == "cruise"
            @assert (prev_stage == "climb" || prev_stage == "loiter" || prev_stage == "descend") "Cruise must be preceded by climb, loiter or descend stage!"
        elseif stage == "loiter"
            @assert (prev_stage == "climb" || prev_stage == "cruise" || prev_stage == "descend") "Loiter must be preceded by climb, cruise or descend stage!"
        elseif stage == "descend"
            @assert (prev_stage == "climb" || prev_stage == "cruise" || prev_stage == "loiter") "Descend must be preceded by climb, cruise or loiter stage!"
            @assert (parse(Float64,df[gradient_row,col]) < 0) "Descend gradient must be smaller than zero!"
        elseif stage == "landing"
            @assert (prev_stage == "descend") "Landing must be preceded by a descend stage!"
        else
            throw(ArgumentError("Invalid Argument: $stage is not a valid stage!"))
        end
    end
end

function convert_to_unitful(; df::DataFrame, N_config::Int)
    # Checks to see whether the inputs are sensisble
    if N_config < 1
        throw(ArgumentError("Number of configurations specified must be bigger than 1!"))
    end

    # For each column
    for col in (ncol(df)-N_config+1):ncol(df)
        # Convert the column to accept anything
        df[!, col] = Vector{Any}(df[!, col])
    end

    # For each input column
    for row in 1:nrow(df)
        # Get the variable type
        type = df[row,"Variable Type"]

        # If the type is not string, then we need to convert
        if type != "String"
            # Parse the variable type
            T = eval(Meta.parse(type))

            # Get the unit
            unit = df[row,"Units"]

            # For each column
            for col in (ncol(df)-N_config+1):ncol(df)
                # If the data is not missing
                if !ismissing(df[row,col])
                    # Convert the data into the given variable type
                    df[row,col] = tryparse(T,string(df[row,col]))

                    # If the unit is defined
                    if !ismissing(unit)
                        # Convert to units for Unitful
                        unit_convert = AeroUnits.convert_to_unit(unit)
                    
                        # Add the units
                        df[row,col] = df[row,col] * unit_convert
                    end
                end
            end
        end
    end

    return df
end

"""
    `validate_data` - A function which validates mission, payload and aircraft user inputs.

    return the cleaned dataframe or throws error
"""
function validate_data(; df :: DataFrame, design_param :: DataFrame, data_type :: String)
    # Checks whether the number of columns make sense
    if ncol(df) <= 2
        throw("Insufficient Data, Please Provide at least 1 Baseline $data_type")
    end

    allowed_type = ["aircraft", "payload", "mission"]

    if !(lowercase(data_type) in allowed_type)
        error("Invalid data_type: $data_type. Must be one of: $(join(allowed_type, ", "))")
    end
    
    # Check that the row and column inputs were unique
    rowcol_names_unique(df)

    # Get the sample data to compare inputs
    df_assumptions = ParseData.get_data(data = "Assumptions");
    df_assumptions = df_assumptions[findall(==(data_type),df_assumptions[!,"Category"]),:]

    # Columns which must exist in an input
    required_cols = ["Assumptions", "Units"];
    df = rowcolumn_check(df = df, df_assumptions = df_assumptions, array_check = required_cols, axis = "Column", mode = "Mandatory")
    
    # To get the number of configurations specified (e.g., aircraft)
    # This is the same as subtracting the number of required columns
    N_config = ncol(df) - length(required_cols)

    # Cheks the row inputs, and append / throws errors if it is not satisfied
    if nrow(df) == 0
        # If the dataframe is empty, then use sample
        @warn "No values were provided, using sample data instead..."
        df = ParseData.get_data(data = "Sample$data_type");
    else
        if data_type == "Mission"
            df = rowcolumn_check(df = df, df_assumptions = df_assumptions, array_check = df_assumptions[:,1], axis = "Row", mode = "Mandatory")

        else
            # If not empty, then check that the row are complete (use df_assumptions to check) and append default values if missing
            df = rowcolumn_check(df = df, df_assumptions = df_assumptions, array_check = df_assumptions[:,1], axis = "Row", mode = "Append", append_rowcol = "Append", verbose = true)
        end
    end

    # List all the columns needed for the validation, which can be appended
    full_cols = vcat(required_cols, ["Category", "Logic", "Mandatory", "Variable Type"])
    # Append the needed columns from the list of assumptions. This SHOULD not require append_rowcol, as all required rows should have been added
    # But it will throw an error if it could not find the assumption
    df = rowcolumn_check(df = df, df_assumptions = df_assumptions, array_check = full_cols, axis = "Column", mode = "Append", append_rowcol = "__", verbose = false)

    # Get the remaining columns which are configurations input (excluding those already in `full_cols`)
    remaining = setdiff(names(df), full_cols)

    # Concatenate to get a final order, with the full_cols at the front and the configurations at the back
    final_order = vcat(full_cols, remaining)

    # Rearrange the column orders
    df = df[:, final_order]

    # Replace names in value inputs into the corresponding input from the referenced column
    if data_type == "Mission"
        df = replace_reference(df = df, key_cols = full_cols, data_type = data_type, skip_row = [1], verbose = false)
    else
        df = replace_reference(df = df, key_cols = full_cols, data_type = data_type, verbose = false)
    end

    # Save any design parameters
    (df, design_param) = design_check(df = df,design_param = design_param, location = data_type, N_config = N_config)

    # Check whether mandatory inputs were specified
    mandatory_check(df = df,key_col = "Mandatory",N_config = N_config)

    # Checks the values of the inputs using the logic column
    df = parse_validation(df = df,key_col = "Logic",N_config = N_config)

    if data_type == "Mission"
        mission_check(df = df, N_config = N_config)
    end

    df = convert_to_unitful(df = df, N_config = N_config)

    # Return the results from the validation and cleaning
    return (df, design_param, N_config)
end

"""
    `df_update_or_append` - A function which update or append dataframe rows with values

    Returns the updated dataframes
"""
function df_update_or_append(;df::DataFrame,label::String,value,N_config::Int,col::Int)
    idx = findfirst(==(label),df[:,1])

    if idx === nothing
        row = deepcopy(df[1, :])

        # Initialise the values as the following
        row["Assumptions"] = label
        row["Units"] = missing
        row["Category"] = "Calculated"
        row["Logic"] = missing
        row["Mandatory"] = true
        row["Variable Type"] = string(eltype(value))

        # Set all the configuration columns as missing
        for idx in (ncol(df)-N_config+1):ncol(df)
            row[idx] = missing
        end

        # Add the initialised column value
        row[col] = value

        # Append the row
        push!(df, row, promote = true)
    else
        # Otherwise, just update/overwrite the value
        df[idx,col] = value
    end

    return df
end

"""
    `update_df_with_design` - A function which updates dataframe with design variables

    Returns the updated dataframes
"""
function update_df_with_design(;design_list::DataFrame,parameter::String,value,df_aircraft::DataFrame,df_mission::DataFrame,df_payload::DataFrame)
    # Find the row corresponding to where the design was saved
    row_idx = findfirst(x -> occursin(Regex("(?i)"*parameter), x), design_list[:,"Design Parameter"])
    location = design_list[row_idx, "Saved Location"]
    row = design_list[row_idx, "Saved Row"]
    col = design_list[row_idx, "Saved Column"]
    
    if location == "Aircraft"
        df_aircraft[row,col] = typeof(df_aircraft[row,col])(value)
    elseif location == "Mission"
        df_mission[row,col] = typeof(df_mission[row,col])(value)
    elseif location == "Payload"
        df_payload[row,col] = typeof(df_payload[row,col])(value)
    else
        throw(ErrorException("Invalid location, cannot find $location"))
    end

    return (df_aircraft,df_mission,df_payload)
end

"""
    `get_design_value` - A function which obtains the design value

    Returns the updated dataframes
"""
function get_design_values(design_list::DataFrame,parameter::String,df_aircraft::DataFrame,df_mission::DataFrame,df_payload::DataFrame)
    # Find the row corresponding to where the design was saved
    row_idx = findfirst(x -> occursin(Regex("(?i)"*parameter), x), design_list[:,"Design Parameter"])
    location = design_list[row_idx, "Saved Location"]
    row = design_list[row_idx, "Saved Row"]
    col = design_list[row_idx, "Saved Column"]
    
    if location == "Aircraft"
        value = df_aircraft[row,col]
    elseif location == "Mission"
        value = df_mission[row,col]
    elseif location == "Payload"
        value = df_payload[row,col]
    else
        throw(ErrorException("Invalid location, cannot find $location"))
    end

    return value
end

end

