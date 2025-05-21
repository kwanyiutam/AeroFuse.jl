module AeroUnits

using Unitful

# Define knots
@unit knots "knots" Knots (1 / 1.94384) * u"m/s" false
# Define horsepower
@unit hp "hp" Horsepower (745.699872) * u"W" false
# Define gallon
@unit gallon "gallon" Gallon (0.0037854118) * u"m^3" false

function convert_to_unit(unit)
    #try
    #    unit = uparse(unit)
    #catch
    #    unit = uparse(unit, unit_context = AeroUnits)
    #end
    # Special case for degrees
    if unit == "degrees"
        unit = "°"
    end

    unit = uparse(unit, unit_context = [Unitful, AeroUnits])
    return unit
end

end
