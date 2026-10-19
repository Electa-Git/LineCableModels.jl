"""
    LineCableModelsMeasurementsExt

Materialize `UncertainValue` as Measurements values and preserve those values
through numerical kernels, display, and data exchange.
"""
module LineCableModelsMeasurementsExt

import Measurements
import Serialization
import SpecialFunctions
using LinearAlgebra: svd
#! explicit-imports: off
# Non-exported accessors documented in Measurements usage.
# Measurements.result is its documented internal propagation
# hook, required by the existing complex-Bessel adapter (Measurements appendix).
# Keep the import exception restricted to these upstream bindings.
using Measurements: value as measured_value, uncertainty as measured_uncertainty,
                    result as measured_result, derivative, uncertainty_components
#! explicit-imports: on

import LineCableModels
import LineCableModels.ParametricBuilder
import LineCableModels.Engine
import LineCableModels.ImportExport

import LineCableModels: nominal, uncertainty
import LineCableModels.Engine: has_uncertainty_type
import LineCableModels.ImportExport:
                                     serialize_value, deserialize_extension,
                                     deserialize_value

# Numeric presentation hooks.
nominal(value::Measurements.Measurement) = measured_value(value)
uncertainty(value::Measurements.Measurement) = measured_uncertainty(value)

# Refinement must still see an uncertain contribution whose nominal value is
# zero. Physical evaluations retain their correlated derivative information.
# Only the scalar numerical error metric and sampling coordinates are nominal.
function Engine.tolerance_norm(z::Complex{<:Measurements.Measurement})
    return max(abs(measured_value(z)),
        hypot(measured_uncertainty(real(z)), measured_uncertainty(imag(z))))
end
function Engine.tolerance_norm(z::Measurements.Measurement)
    return max(abs(measured_value(z)), measured_uncertainty(z))
end

function LineCableModels.materialize(value::ParametricBuilder.UncertainValue{<:Real})
    Measurements.measurement(value.nominal, value.sigma)
end

function has_uncertainty_type(
        ::Type{Complex{T}},
) where {T <: Measurements.Measurement}
    true
end
has_uncertainty_type(::Type{<:Measurements.Measurement}) = true
function serialize_value(value::Measurements.Measurement)
    return Dict(
        "__type__" => "Measurement",
        "value" => serialize_value(measured_value(value)),
        "uncertainty" => serialize_value(measured_uncertainty(value))
    )
end
function deserialize_extension(::Val{:Measurement}, value)
    nominal = deserialize_value(value["value"])
    uncertainty = deserialize_value(value["uncertainty"])
    return Measurements.measurement(nominal, uncertainty)
end

# A loaded result gives its uncertainty sources new identities, one per source of the file.
# Its values keep their correlations with one another and are independent of this process's
# sources. A dependent value keeps its nominal value and uncertainty, and its derivatives keep
# their order. `measured_result` rebuilds the derivatives against the new sources. `tag` and
# `der` are the fields that the Measurements appendix documents.
function Serialization.deserialize(reader::ImportExport.ResultReader,
        type::Type{<:Measurements.Measurement})
    saved = invoke(Serialization.deserialize,
        Tuple{Serialization.AbstractSerializer, DataType}, reader, type)
    source(key) = get!(() -> Measurements.measurement(key[1], key[2]), reader.sources, key)
    saved.tag == 0 || return source((measured_value(saved), measured_uncertainty(saved), saved.tag))
    # `der` iterates from its newest entry.
    derivatives = reverse!(collect(saved.der))
    isempty(derivatives) && return saved
    rebuilt = measured_result(measured_value(saved), last.(derivatives), source.(first.(derivatives)))
    return Measurements.Measurement(measured_value(saved), measured_uncertainty(saved), UInt64(0),
        rebuilt.der)
end
# Uncertainty-aware SpecialFunctions methods used by the numerical kernels.
function _lift_complex(function_value, order, value::Complex{<:Measurements.Measurement})
    z = measured_value(value)
    result = function_value(order, z)
    lower, upper = function_value(order - 1, z), function_value(order + 1, z)
    # DLMF 10.6.1 and 10.29.1: exact order recurrences avoid an absolute
    # finite-difference step crossing the singularity at small arguments.
    derivative = if function_value in (SpecialFunctions.besseli, SpecialFunctions.besselix)
        (lower + upper) / 2
    elseif function_value in (SpecialFunctions.besselk, SpecialFunctions.besselkx)
        -(lower + upper) / 2
    else
        (lower - upper) / 2
    end
    dx, dy = derivative, im * derivative
    # Scaled I/J/Y are not holomorphic. Differentiate their scaling in Cartesian
    # coordinates. `sign(0)=0` retains the symmetric slope at an absolute-value cusp.
    if function_value === SpecialFunctions.besselix
        dx -= sign(real(z)) * result
    elseif function_value in (SpecialFunctions.besseljx, SpecialFunctions.besselyx)
        dy -= sign(imag(z)) * result
    elseif function_value === SpecialFunctions.besselkx
        dx += result
        dy = im * dx
    elseif function_value === SpecialFunctions.besselhx
        dx -= im * result
        dy = im * dx
    end
    return measured_result(result, [real(dx), real(dy), imag(dx), imag(dy)], value)
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

"""
Lift local capacitance through its fixed-resolution least-squares equations.
Independent uncertainty directions preserve the original Measurement graph.
Only the nominal factorization and bounded Float64 derivative blocks are dense.
Kernel derivatives use centered steps. This avoids both a dense Measurement matrix
and repeated perturbed factorizations.
"""
function Engine.ShuntModel.internal_capacitance(selected::Engine.ShuntModel.Formula{:boundary},
        group::Engine.ShuntModel.DielectricGroup{<:Measurements.Measurement})
    values = Vector(group)
    nominal_group = Engine.ShuntModel.DielectricGroup(group, Float64.(measured_value.(values)))
    tags = sort!(
        unique!([tag for value in values
                 for tag in keys(uncertainty_components(value)) if !iszero(tag[2])]);
        by = last)
    if isempty(tags)
        result = Engine.ShuntModel.internal_capacitance(selected, nominal_group)
        return (; C = eltype(values).(result.C), result.diagnostic)
    end
    # Columns are physical perturbations per standard deviation of each
    # independent input. Their distribution remains unchanged by this normalization.
    directions = [Float64(derivative(value, tag)*tag[2]) for value in values, tag in tags]
    result = Engine.ShuntModel.internal_capacitance(selected, nominal_group, directions)
    # Return sensitivities through the original physical arguments, including
    # shared and dependent ones. A small rank-revealing solve handles redundant
    # descriptors without creating new independent Measurement identities.
    factor = svd(transpose(directions))
    cutoff = max(size(directions)...)*eps(Float64)*maximum(factor.S)
    keep = factor.S .> cutoff
    gradient = factor.V[:, keep]*((transpose(factor.U[:, keep])*result.tangents) ./
                                  factor.S[keep])
    lifted = [measured_result(result.C[i], @view(gradient[:, i]), values)
              for i in eachindex(result.C)]
    return (; C = reshape(lifted, size(result.C)), result.diagnostic)
end
# The archive source table spans all atomic points and a separate report reference.
function ImportExport.encode_observation(value::Measurements.Measurement,source_table)
    contributions=NamedTuple[]
    for source in keys(uncertainty_components(value))
        index=get!(source_table.indices,source) do
            push!(source_table.sources,(nominal=source[1],sigma=source[2]))
            length(source_table.sources)
        end
        push!(contributions,(source=index,sensitivity=derivative(value,source)))
    end
    sort!(contributions;by=entry -> entry.source)
    return Dict("__type__"=>"ObservedMeasurement",
        "nominal"=>serialize_value(nominal(value)),
        "uncertainty"=>serialize_value(uncertainty(value)),
        "contributions"=>serialize_value(contributions,Val(:scientific)))
end
function ImportExport.observation_sources(records::AbstractVector,::Val{:measurements})
    return [Measurements.measurement(record.nominal,record.sigma) for record in records]
end
function ImportExport.decode_observation_measurement(record,sources)
    center=deserialize_value(record["nominal"])
    value=Measurements.measurement(center,zero(center))
    for entry in deserialize_value(record["contributions"])
        source=sources[entry.source]
        value+=entry.sensitivity*(source-nominal(source))
    end
    expected=deserialize_value(record["uncertainty"])
    isapprox(uncertainty(value),expected;rtol=1e-12,atol=0) ||
        throw(ArgumentError("saved uncertainty sensitivities do not reproduce the retained standard deviation"))
    return value
end

end
