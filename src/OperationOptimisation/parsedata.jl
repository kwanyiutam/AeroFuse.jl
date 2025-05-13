module ParseData

# Initialise packages used
using CSV
using DataFrames

"""
    `get_data` - A function which obtains defined default tables from CSV files.

    Takes `data::String`, which corresponds to a dataset ID, and the function will try and find the dataset under the main csv.

    return the dataframe of the CSV file
"""
function get_data(; data :: String)
    # Get the master dataset
    master = CSV.read("MasterData.csv", DataFrame; stringtype=String)

    # Try and find the dataset with the matching "data" ID 
    result = try
        master[only(findall(==(data), master.Data)), "Path"]
    catch e
        throw("$data dataset cannot be found!")
    end

    # Get the dataframe result
    df_result = CSV.read(result, DataFrame; stringtype=String)

    # Return the results from the search
    return df_result
end

end