module LineCableModelsMeasurementsExt

using Calculus
using Measurements
using Printf
using SpecialFunctions

import LineCableModels
const Grammar = LineCableModels.Grammar
const ParametricBuilder = LineCableModels.ParametricBuilder
const Engine = LineCableModels.Engine
const ImportExport = LineCableModels.ImportExport

import LineCableModels: nominal, standard_uncertainty

# Numeric presentation hooks.
nominal(value::Measurements.Measurement) = Measurements.value(value)
standard_uncertainty(value::Measurements.Measurement) = Measurements.uncertainty(value)

function ParametricBuilder.materialize(value::ParametricBuilder.UncertainValue{<:Real})
    Measurements.measurement(value.nominal, value.sigma)
end

function Engine._has_uncertainty_type(
        ::Type{Complex{T}},
) where {T <: Measurements.Measurement}
    true
end
function Engine._clip_field(value::Measurements.Measurement, tolerance)
    nominal = abs(Measurements.value(value)) <= tolerance ?
              0.0 : Measurements.value(value)
    uncertainty = abs(Measurements.uncertainty(value)) <= tolerance ?
                  0.0 : Measurements.uncertainty(value)
    return Measurements.measurement(nominal, uncertainty)
end

function ImportExport._serialize_value(value::Measurements.Measurement)
    return Dict(
        "__type__" => "Measurement",
        "value" => ImportExport._serialize_value(Measurements.value(value)),
        "uncertainty" => ImportExport._serialize_value(Measurements.uncertainty(value))
    )
end
function ImportExport._deserialize_extension(::Val{:Measurement}, value)
    nominal = ImportExport._deserialize_value(value["value"])
    uncertainty = ImportExport._deserialize_value(value["uncertainty"])
    return Measurements.measurement(nominal, uncertainty)
end
function ImportExport.stringify(value::Measurements.Measurement)
    Printf.@sprintf("%.12g ± %.6g",
        Measurements.value(value),
        Measurements.uncertainty(value),)
end

# Uncertainty-aware SpecialFunctions methods used by the numerical kernels.
function _lift_complex(function_value, order, value::Complex{<:Measurements.Measurement})
    nominal = Measurements.value(value)
    return Measurements.result(
        function_value(order, nominal),
        vcat(
            Calculus.gradient(
                point -> real(function_value(order, complex(point...))),
                collect(reim(nominal))
            ),
            Calculus.gradient(
                point -> imag(function_value(order, complex(point...))),
                collect(reim(nominal))
            )
        ),
        value
    )
end

for name in (
    :besselix, :besselkx, :besseljx, :besselyx, :besselhx,
    :besseli, :besselk, :besselj, :bessely, :besselh
)
    @eval begin
        function SpecialFunctions.$name(
                order::Real,
                value::Complex{<:Measurements.Measurement}
        )
            return _lift_complex(SpecialFunctions.$name, order, value)
        end
    end
end

end
