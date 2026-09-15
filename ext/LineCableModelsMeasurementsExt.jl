"""
    LineCableModelsMeasurementsExt

Materialize `UncertainValue` as Measurements values and preserve those values
through numerical kernels, display, and data exchange.
"""
module LineCableModelsMeasurementsExt

import Measurements
import Printf
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
import LineCableModels.UQ
import LineCableModels.ReportBuilder

import LineCableModels: nominal, uncertainty
import LineCableModels.Engine: has_uncertainty_type
import LineCableModels.ImportExport:
                                     serialize_value, deserialize_extension,
                                     deserialize_value
import LineCableModels.Grammar: detach
import LineCableModels.ReportBuilder: encode_cell

# Numeric presentation hooks.
nominal(value::Measurements.Measurement) = measured_value(value)
uncertainty(value::Measurements.Measurement) = measured_uncertainty(value)

# Refinement must still see an uncertain contribution whose nominal value is
# zero. Physical evaluations retain their correlated derivative information;
# only the scalar numerical error metric and sampling coordinates are nominal.
function Engine.spectral_magnitude(z::Complex{<:Measurements.Measurement})
    return max(abs(measured_value(z)),
        hypot(measured_uncertainty(real(z)), measured_uncertainty(imag(z))))
end
function Engine.spectral_magnitude(z::Measurements.Measurement)
    return max(abs(measured_value(z)), measured_uncertainty(z))
end

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
        representative.frequency,
        representative.details
    )
end

function _measurement_result(
        source::UQ.MonteCarloResult{<:Engine.LineParameters{T, U, D, Basis}},
        point::Integer
) where {T, U, D, Basis}
    representative = source.values[point]
    summary = source.stats[point]
    resistance = _measurement.(summary.R)
    inductance = _measurement.(summary.L)
    capacitance = _measurement.(summary.C)
    conductance = _measurement.(summary.G)
    angular = reshape(2π .* representative.f, 1, 1, :)
    impedance = complex.(resistance, inductance .* angular)
    admittance = complex.(conductance, capacitance .* angular)
    element_type = promote_type(eltype(impedance), eltype(admittance))
    return Engine.LineParameters(
        representative.domain,
        Engine.SeriesImpedance{element_type, Basis}(convert(Array{element_type, 3}, impedance)),
        Engine.ShuntAdmittance{element_type, Basis}(convert(Array{element_type, 3}, admittance)),
        representative.f,
        representative.details
    )
end

function UQ.uncertain(source::UQ.MonteCarloResult{T}, point::Integer) where {
        T <: Union{Engine.CableConstants, Engine.LineParameters}}
    checkbounds(source.values, point)
    return _measurement_result(source, point)
end

function UQ.uncertain(
        source::UQ.MonteCarloResult{T}
) where {
        T <: Union{Engine.CableConstants, Engine.LineParameters}
}
    first_value = UQ.uncertain(source, firstindex(source))
    values = Vector{typeof(first_value)}(undef, length(source))
    values[1] = first_value
    for point in 2:length(source)
        value = UQ.uncertain(source, point)
        typeof(value) === eltype(values) || throw(ArgumentError(
            "Monte Carlo statistics reconstructed inconsistent result types",
        ))
        values[point] = value
    end
    return values
end

function ParametricBuilder.Gridspace{Target}(source::UQ.MonteCarloResult{T}) where {
        Target, T <: Union{Engine.CableConstants, Engine.LineParameters}}
    return ParametricBuilder.Gridspace{Target}(Target, (ParametricBuilder.Grid(UQ.uncertain(source)),))
end

function has_uncertainty_type(
        ::Type{Complex{T}},
) where {T <: Measurements.Measurement}
    true
end
function detach(value::Measurements.Measurement, factor, clip::Bool)
    return value * factor
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
        "value" => serialize_value(measured_value(value)),
        "uncertainty" => serialize_value(measured_uncertainty(value))
    )
end
function deserialize_extension(::Val{:Measurement}, value)
    nominal = deserialize_value(value["value"])
    uncertainty = deserialize_value(value["uncertainty"])
    return Measurements.measurement(nominal, uncertainty)
end

"""
Encode one LEP result with shared independent-source identities and sparse signed
sensitivities. Source IDs are local to this record; separate outputs are never
reconstructed as independent measurements.
"""
function serialize_value(value::UQ.LinearErrorResult{<:Engine.LineParameters{<:Complex{<:Measurements.Measurement}}})
    record=NamedTuple(value)
    isempty(value) && throw(ArgumentError("cannot encode an empty LEP result"))
    initial=real(first(LineCableModels.observe(first(value),LineCableModels.Z)))
    source_set=Set(keys(uncertainty_components(initial)))
    for point_result in value, selector in (LineCableModels.Z,LineCableModels.Y),
            z in LineCableModels.observe(point_result,selector), component in (real(z),imag(z))
        union!(source_set,keys(uncertainty_components(component)))
    end
    sources=sort!(collect(source_set);by=last)
    indices=Dict(source=>index for (index,source) in enumerate(sources))
    points=map(value) do core
        matrices=map((LineCableModels.Z,LineCableModels.Y)) do selector
            values=LineCableModels.observe(core,selector)
            components=map((real,imag)) do component_part
                [begin
                    component=component_part(z)
                    contributions=[(source=indices[source],sensitivity=derivative(component,source))
                        for source in keys(uncertainty_components(component))]
                    sort!(contributions;by=entry -> entry.source)
                    (nominal=nominal(component),uncertainty=uncertainty(component),contributions)
                end for z in values]
            end
            (shape=size(values),real=vec(components[1]),imaginary=vec(components[2]))
        end
        (Z=matrices[1],Y=matrices[2],frequencies=LineCableModels.frequencies(core),
            basis=LineCableModels.basis(core),domain=:PhaseDomain,
            coordinates=get(LineCableModels.details(core),:coordinates,nothing),
            shunt_model=get(LineCableModels.details(core),:shunt_model,nothing))
    end
    formulation=record.formulation isa NamedTuple ? record.formulation : NamedTuple(record.formulation)
    payload=(formulation,points,sources=[(nominal=source[1],sigma=source[2]) for source in sources],details=record.details)
    return Dict("__type__"=>"MeasurementLinearErrorResult","version"=>1,
        "payload"=>serialize_value(payload,Val(:scientific)))
end

function serialize_value(value::UQ.LinearErrorResult{<:Engine.CableConstants{<:Measurements.Measurement}})
    isempty(value) && throw(ArgumentError("cannot encode an empty LEP result"))
    initial=first(first(value).R)
    source_set=Set(keys(uncertainty_components(initial)))
    for core in value, component in Iterators.flatten((core.R,core.L,core.C,core.G,(core.frequency,)))
        union!(source_set,keys(uncertainty_components(component)))
    end
    sources=sort!(collect(source_set);by=last)
    indices=Dict(source=>index for (index,source) in enumerate(sources))
    encode = component -> begin
        contributions=[(source=indices[source],sensitivity=derivative(component,source))
            for source in keys(uncertainty_components(component))]
        sort!(contributions;by=entry->entry.source)
        (nominal=nominal(component),uncertainty=uncertainty(component),contributions)
    end
    points=map(value) do core
        (kind=:cable_constants,cores=core.cores,R=encode.(core.R),L=encode.(core.L),
            C=encode.(core.C),G=encode.(core.G),frequency=encode(core.frequency),
            shunt_model=get(LineCableModels.details(core),:shunt_model,nothing))
    end
    record=NamedTuple(value)
    formulation=record.formulation isa NamedTuple ? record.formulation : NamedTuple(record.formulation)
    payload=(formulation,points,sources=[(nominal=source[1],sigma=source[2]) for source in sources],details=record.details)
    return Dict("__type__"=>"MeasurementLinearErrorResult","version"=>1,
        "payload"=>serialize_value(payload,Val(:scientific)))
end

function deserialize_extension(::Val{:MeasurementLinearErrorResult},record)
    record["version"] == 1 || throw(ArgumentError("unsupported shared-source LEP record"))
    payload=deserialize_value(record["payload"])
    sources=[Measurements.measurement(source.nominal,source.sigma) for source in payload.sources]
    restore = component -> begin
        value=Measurements.measurement(component.nominal,zero(component.nominal))
        for entry in component.contributions
            source=sources[entry.source]
            value += entry.sensitivity*(source-nominal(source))
        end
        isapprox(uncertainty(value),component.uncertainty;rtol=1e-12,atol=0) ||
            throw(ArgumentError("restored LEP sensitivity record changes propagated uncertainty"))
        value
    end
    points=map(payload.points) do point
        model=get(point,:shunt_model,nothing)
        detail=model===nothing ? (;) : (shunt_model=model,)
        if get(point,:kind,nothing)===:cable_constants
            return Engine.CableConstants(point.cores,restore.(point.R),restore.(point.L),
                restore.(point.C),restore.(point.G),restore(point.frequency),detail)
        end
        matrices=map((point.Z,point.Y)) do matrix
            parts=map((matrix.real,matrix.imaginary)) do components
                restore.(components)
            end
            reshape(complex.(parts...),matrix.shape)
        end
        point.domain === :PhaseDomain || throw(ArgumentError("unsupported LEP result domain"))
        point.coordinates === nothing || (detail=merge(detail,(coordinates=point.coordinates,)))
        Engine.LineParameters(matrices...,point.frequencies;basis=point.basis,details=detail)
    end
    return UQ.LinearErrorResult(payload.formulation,points,payload.details)
end
function encode_cell(
        ::ReportBuilder.XLSXReportDefinition,
        value::Measurements.Measurement
)
    Printf.@sprintf("%.12g ± %.6g",
        measured_value(value),
        measured_uncertainty(value),)
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
    # coordinates; sign(0)=0 retains the symmetric slope at an absolute-value cusp.
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
Independent uncertainty directions preserve the original Measurement graph;
only the nominal factorization and bounded Float64 derivative blocks are dense.
Kernel derivatives use centered steps, checked by halving the step. This avoids
both a dense Measurement matrix and repeated perturbed factorizations.
"""
function Engine.internal_shunt_response(values::AbstractVector{<:Measurements.Measurement},domain;kwargs...)
    nominal_values = Float64.(measured_value.(values))
    tags = sort!(unique!([tag for value in values
        for tag in keys(uncertainty_components(value)) if !iszero(tag[2])]);by=last)
    if isempty(tags)
        result = Engine.internal_shunt_response(nominal_values,domain;kwargs...)
        return (;C=eltype(values).(result.C),diagnostic=result.diagnostic,state=nothing)
    end
    # Columns are physical perturbations per standard deviation of each
    # independent input. This normalization does not change their distribution.
    directions = [Float64(derivative(value,tag)*tag[2]) for value in values, tag in tags]
    result = Engine.internal_shunt_response(nominal_values,domain;directions,kwargs...)
    # Return sensitivities through the original physical arguments, including
    # shared/dependent ones. A small rank-revealing solve handles redundant
    # descriptors without creating new independent Measurement identities.
    factor = svd(transpose(directions))
    cutoff = max(size(directions)...)*eps(Float64)*maximum(factor.S)
    keep = factor.S .> cutoff
    gradient = factor.V[:,keep]*((transpose(factor.U[:,keep])*result.tangents)./factor.S[keep])
    lifted = [measured_result(result.C[i],@view(gradient[:,i]),values) for i in eachindex(result.C)]
    return (;C=reshape(lifted,size(result.C)),diagnostic=result.diagnostic,state=nothing)
end
end
