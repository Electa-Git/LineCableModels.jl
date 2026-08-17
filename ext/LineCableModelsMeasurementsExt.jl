module LineCableModelsMeasurementsExt

using Calculus
using Measurements
using Printf
using SpecialFunctions
using Statistics

import LineCableModels
const PB = LineCableModels.ParametricBuilder
const Engine = LineCableModels.Engine
const DataModel = LineCableModels.DataModel
const ImportExport = LineCableModels.ImportExport
const Computation = LineCableModels.Computation

import LineCableModels: nominal, standard_uncertainty

# Numeric presentation hooks.
nominal(value::Measurements.Measurement) = Measurements.value(value)
standard_uncertainty(value::Measurements.Measurement) = Measurements.uncertainty(value)

# Gridspace direct propagation. The realization cache in Gridspace guarantees
# that one shared Grid becomes one shared Measurement variable, including
# across nested object boundaries.
function PB._direct_value(value::PB.UncertainValue{<:Real})
    Measurements.measurement(value.nominal, value.sigma)
end
PB.materialize(value::PB.UncertainValue{<:Real}) = PB._direct_value(value)

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

function _joint_coordinates(matrix::AbstractMatrix{T}) where {T <: Real}
    means = vec(Statistics.mean(matrix; dims = 2))
    size(matrix, 2) == 1 && return map(
        value -> Measurements.measurement(value, zero(T)), means
    )
    factor = matrix .- means
    factor ./= sqrt(size(matrix, 2) - 1)
    latent = [Measurements.measurement(zero(T), one(T)) for _ in axes(matrix, 2)]
    return means + factor * latent
end

function _independent_coordinates(means, deviations)
    return map(Measurements.measurement, means, deviations)
end

function Measurements.measurement(result::Computation.MonteCarloResult{<:DataModel.CableConstants})
    if result.samples === nothing
        @warn "MonteCarloResult has no retained samples; reconstructing independent marginal Measurements and discarding output correlations"
        mean = result.representation
        sigma = result.uncertain.sigma
        return DataModel.CableConstants(
            Measurements.measurement(mean.R, sigma.R),
            Measurements.measurement(mean.L, sigma.L),
            Measurements.measurement(mean.C, sigma.C)
        )
    end
    matrix = permutedims(hcat(
        result.samples.R,
        result.samples.L,
        result.samples.C
    ))
    return DataModel.CableConstants(_joint_coordinates(matrix)...)
end

function _line_measurements(result::Computation.MonteCarloResult{<:Engine.LineParameters})
    nominal = result.representation
    tensor_size = size(nominal.Z)
    coordinate_count = prod(tensor_size)
    if result.samples === nothing
        @warn "MonteCarloResult has no retained samples; reconstructing independent marginal Measurements and discarding output correlations"
        means = (
            real.(nominal.Z.values),
            imag.(nominal.Z.values) ./ reshape(2π .* nominal.f, 1, 1, :),
            real.(nominal.Y.values),
            imag.(nominal.Y.values) ./ reshape(2π .* nominal.f, 1, 1, :)
        )
        deviations = result.uncertain.sigma
        blocks = map(
            _independent_coordinates,
            means,
            (deviations.R, deviations.L, deviations.G, deviations.C)
        )
    else
        sample_values = result.samples
        trial_count = size(sample_values.R, 4)
        matrix = vcat(
            reshape(sample_values.R, coordinate_count, trial_count),
            reshape(sample_values.L, coordinate_count, trial_count),
            reshape(sample_values.G, coordinate_count, trial_count),
            reshape(sample_values.C, coordinate_count, trial_count)
        )
        joint = _joint_coordinates(matrix)
        blocks = ntuple(
            index -> reshape(
                joint[((index - 1) * coordinate_count + 1):(index * coordinate_count)],
                tensor_size
            ),
            4)
    end
    resistance, inductance, conductance, capacitance = blocks
    angular = reshape(2π .* nominal.f, 1, 1, :)
    impedance = complex.(resistance, inductance .* angular)
    admittance = complex.(conductance, capacitance .* angular)
    return Engine.LineParameters(
        LineCableModels.domain(nominal),
        impedance,
        admittance,
        nominal.f;
        basis = LineCableModels.basis(nominal)
    )
end

function Measurements.measurement(result::Computation.MonteCarloResult{<:Engine.LineParameters})
    _line_measurements(result)
end

end
