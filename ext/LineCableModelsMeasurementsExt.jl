"""
    LineCableModelsMeasurementsExt

Materialise `UncertainValue` as Measurements values and preserve those values
through numerical kernels, display, and data exchange.
"""
module LineCableModelsMeasurementsExt

using Calculus
using Measurements
using Printf
using SpecialFunctions

import LineCableModels
const Grammar = LineCableModels.Grammar
const ParametricBuilder = LineCableModels.ParametricBuilder
const Engine = LineCableModels.Engine
const UQ = LineCableModels.UQ
const ReportBuilder = LineCableModels.ReportBuilder
const ImportExport = LineCableModels.ImportExport

import LineCableModels: nominal, uncertainty
import LineCableModels.Engine: has_uncertainty_type
import LineCableModels.ImportExport:
                                     serialize_value, deserialize_extension,
                                     deserialize_value
import LineCableModels.Grammar: detach
import LineCableModels.ReportBuilder: encode_cell

# Numeric presentation hooks.
nominal(value::Measurements.Measurement) = Measurements.value(value)
uncertainty(value::Measurements.Measurement) = Measurements.uncertainty(value)

function LineCableModels.materialize(value::ParametricBuilder.UncertainValue{<:Real})
    Measurements.measurement(value.nominal, value.sigma)
end

function _measurement(summary::UQ.SampleSummary)
    Measurements.measurement(summary.mean, summary.std)
end

function _measurement_result(
        source::UQ.MonteCarloResult{<:Engine.CableConstants},
        point::Integer
)
    representative = source.values[point]
    summary = source.stats[point]
    return Engine.CableConstants(
        representative.cores,
        _measurement.(summary.R),
        _measurement.(summary.L),
        _measurement.(summary.C),
        _measurement.(summary.G),
        representative.frequency
    )
end

function _measurement_result(
        source::UQ.MonteCarloResult{<:Engine.LineParameters},
        point::Integer
)
    representative = source.values[point]
    summary = source.stats[point]
    resistance = _measurement.(summary.R)
    inductance = _measurement.(summary.L)
    capacitance = _measurement.(summary.C)
    conductance = _measurement.(summary.G)
    angular = reshape(2π .* representative.f, 1, 1, :)
    impedance = complex.(resistance, inductance .* angular)
    admittance = complex.(conductance, capacitance .* angular)
    return Engine.LineParameters(
        representative.domain,
        impedance,
        admittance,
        representative.f;
        basis = LineCableModels.basis(representative),
        details = representative.details
    )
end

function ParametricBuilder.Gridspace{Target}(
        source::UQ.MonteCarloResult{T}
) where {
        Target,
        T <: Union{Engine.CableConstants, Engine.LineParameters}
}
    first_value = _measurement_result(source, firstindex(source))
    values = Vector{typeof(first_value)}(undef, length(source))
    values[1] = first_value
    for point in 2:length(source)
        value = _measurement_result(source, point)
        typeof(value) === eltype(values) || throw(ArgumentError(
            "Monte Carlo statistics reconstructed inconsistent result types",
        ))
        values[point] = value
    end
    return ParametricBuilder.Gridspace{Target}(
        Target,
        (ParametricBuilder.Grid(values),)
    )
end

function has_uncertainty_type(
        ::Type{Complex{T}},
) where {T <: Measurements.Measurement}
    true
end
function detach(value::Measurements.Measurement, factor, clip::Bool)
    nominal = detach(Measurements.value(value), factor, clip)
    uncertainty = detach(Measurements.uncertainty(value), abs(factor), clip)
    return Measurements.measurement(nominal, uncertainty)
end

function detach(
        values::AbstractArray{<:Measurements.Measurement},
        factor,
        clip::Bool
)
    return map(value -> detach(value, factor, clip), values)
end

function serialize_value(value::Measurements.Measurement)
    return Dict(
        "__type__" => "Measurement",
        "value" => serialize_value(Measurements.value(value)),
        "uncertainty" => serialize_value(Measurements.uncertainty(value))
    )
end
function deserialize_extension(::Val{:Measurement}, value)
    nominal = deserialize_value(value["value"])
    uncertainty = deserialize_value(value["uncertainty"])
    return Measurements.measurement(nominal, uncertainty)
end
function encode_cell(
        ::ReportBuilder.XLSXReportDefinition,
        value::Measurements.Measurement
)
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

function SpecialFunctions.besselix(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besselix, order, value)
end

function SpecialFunctions.besselkx(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besselkx, order, value)
end

function SpecialFunctions.besseljx(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besseljx, order, value)
end

function SpecialFunctions.besselyx(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besselyx, order, value)
end

function SpecialFunctions.besselhx(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besselhx, order, value)
end

function SpecialFunctions.besseli(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besseli, order, value)
end

function SpecialFunctions.besselk(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besselk, order, value)
end

function SpecialFunctions.besselj(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besselj, order, value)
end

function SpecialFunctions.bessely(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.bessely, order, value)
end

function SpecialFunctions.besselh(order::Real, value::Complex{<:Measurements.Measurement})
    return _lift_complex(SpecialFunctions.besselh, order, value)
end

end
