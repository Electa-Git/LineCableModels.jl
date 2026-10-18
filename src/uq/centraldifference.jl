"""
Store a failed central perturbation and its input diagnostics. Access the
input path, signed step, computation stage and original error with `details`.
"""
struct PerturbationError{D,E} <: Exception
    "Input and perturbation diagnostics."
    details::D
    "Original construction, computation or mode-assignment failure."
    cause::E
end

details(error::PerturbationError) = error.details
function Base.showerror(io::IO, error::PerturbationError)
    record = details(error).data
    print(io, "LinearError central difference failed for ", record.input,
        " at signed step ", record.step, " during ", record.stage, ": ")
    showerror(io, error.cause)
end

function _central_failure(exception::InterruptException, diagnostic)
    rethrow()
end
function _central_failure(exception, diagnostic)
    throw(PerturbationError(ComputationDetails(diagnostic), exception))
end

# Descriptors at Gridpoint leaves are independent. Existing uncertain numerical
# objects carry a source graph that this declaration cannot represent.
function validate(value::Number, ::CentralDifferences)
    Engine.has_uncertainty_type(typeof(value)) && throw(ArgumentError(
        "LinearError central differences cannot represent supplied shared or correlated " *
        "uncertain inputs; declare independent inputs with Grid and derive dependencies in its builder"))
    return value
end
function validate(value::UncertainValue, ::CentralDifferences)
    throw(ArgumentError("LinearError central differences require real independent " *
        "uncertainty descriptors at Gridspace argument leaves; got $(typeof(value))"))
end
function validate(value::Union{AbstractArray,Tuple,NamedTuple,AbstractDict},
        selected::CentralDifferences)
    foreach(entry -> validate(entry, selected), value)
    return value
end
function validate(value::Union{Function,Type,Module}, ::CentralDifferences)
    return value
end
function validate(value, selected::CentralDifferences)
    for name in fieldnames(typeof(value))
        validate(getfield(value, name), selected)
    end
    return value
end

function _central_inputs(method::CentralDifferences, records, point::Gridpoint, path)
    for (index, argument) in enumerate(point.args)
        _central_inputs(method, records, argument, (path..., index))
    end
    return records
end
function _central_inputs(method::CentralDifferences, records, value::UncertainValue{<:Real,<:Real}, path)
    iszero(value.sigma) || push!(records, (input="input" * join("[$index]" for index in path),
        argument_path=path, nominal=value.nominal, standard_uncertainty=value.sigma))
    return records
end
function _central_inputs(method::CentralDifferences, records, value, path)
    validate(value, method)
    return records
end

function _linear_error(method::CentralDifferences, formulation::LinearError, problem,
        modal, modal_options)
    composed = ModalComposition(formulation.inner, modal, modal_options)
    return compute(problem, LinearError(composed, method, formulation.options))
end

# The modal reference is the already computed nominal point, with no extra solve.
function _sample_axis(value::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain},
        point, formulation::LinearError{<:ModalComposition}, child_options)
    selected = formulation.inner.modal.formula
    frequency = Engine.ModalAnalysis.ordering_frequency(selected, frequencies(value))
    return (; reference=value, selected, frequency, frequencies=frequencies(value))
end
function _sample_axis(value::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain},
        point, formulation::LinearError{<:LineCableModels.ModalAnalysisFormulation}, child_options)
    selected = formulation.inner.formula
    frequency = Engine.ModalAnalysis.ordering_frequency(selected, frequencies(value))
    return (; reference=value, selected, frequency, frequencies=frequencies(value))
end

# One independent unit-uncertainty source per input preserves the signed Jacobian
# across all outputs, including real and imaginary components and modal products.
function _central_values(values::AbstractArray, sources, relative_step)
    shape = size(values)[1:end-1]
    return map(CartesianIndices(shape)) do index
        samples = @view values[index.I..., :]
        all(isfinite, samples) || return missing
        center = samples[1]
        result = center + materialize(UncertainValue(zero(real(center)), zero(real(center))))
        for (input, source) in enumerate(sources)
            result += (samples[2input] - samples[2input+1]) / (2relative_step) * source
        end
        result
    end
end
function _central_values(products::NamedTuple, sources, relative_step)
    return map(product -> _central_values(product, sources, relative_step), products)
end

function _central_difference(method::CentralDifferences, point, formulation, child_options, point_index)
    inputs = _central_inputs(method, NamedTuple[], point, ())
    nominal_problem = realize(Random.Xoshiro(0), point, (_rng, center, _sigma) -> float(center))
    reference = compute(nominal_problem, formulation.inner; options=child_options)
    axis = _sample_axis(reference, point, formulation, child_options)
    count = 2length(inputs) + 1
    buffers = initialize_buffers(formulation, Float64, reference, (; trials=count), (;))
    _record_sample!(buffers, reference, 1, axis)
    retained = formulation.options.data.retain_details ?
        [computation_details(typeof(formulation.inner), reference)] : nothing
    relative_step = method.relative_step
    perturbations = NamedTuple[]
    for (input_index, input) in enumerate(inputs)
        step = relative_step * input.standard_uncertainty
        for (side, direction) in enumerate((1, -1))
            stage = :build
            evaluation = 2input_index + side - 1
            try
                isfinite(step) && step > 0 && input.nominal + step != input.nominal &&
                    input.nominal - step != input.nominal || throw(ArgumentError(
                    "the requested perturbation must be finite, positive and distinguishable from the nominal input"))
                leaf = Ref(0)
                shifted = realize(Random.Xoshiro(0), point, (_rng, center, _sigma) -> begin
                    leaf[] += 1
                    float(center) + (leaf[] == input_index ? direction * step : zero(step))
                end)
                stage = :compute
                value = compute(shifted, formulation.inner; options=child_options)
                _sample_shunt_model(value) == _sample_shunt_model(reference) || throw(ArgumentError(
                    "central perturbations changed shunt-model coverage; select a fixed geometry model"))
                stage = :record
                accepted = _record_sample!(buffers, value, evaluation, axis)
                stage = :assignment
                accepted === false && throw(ArgumentError(
                    "modal current-basis overlap $(buffers.overlaps[:,evaluation]) falls below " *
                    "the mode-assignment threshold $_MODAL_OVERLAP_THRESHOLD; " *
                    "central differences require both perturbations"))
                retained === nothing || push!(retained, computation_details(typeof(formulation.inner), value))
                push!(perturbations, (; input..., step=direction*step, evaluation))
            catch exception
                diagnostic = (; point=point_index, input..., step=direction*step, evaluation, stage,
                    mode_assignment=stage === :assignment && haskey(buffers, :overlaps) ?
                        (overlap=copy(buffers.overlaps[:,evaluation]),
                            permutation=copy(buffers.permutations[:,evaluation])) : nothing)
                _central_failure(exception, diagnostic)
            end
        end
    end
    sources = [materialize(UncertainValue(0.0, 1.0)) for _ in inputs]
    products = _central_values(get(buffers, :products, buffers), sources, relative_step)
    representation = materialize(reference, products, formulation)
    representation = Engine.retain_gridpoint(representation, Commons.gridpoint_id(); fields=(
        inputs=merge(get(details(reference).data, :inputs, (;)),
            (interpretation=:nominal_declaration, uncertain_arguments=_uncertain_arguments(point))),
        uncertainty=(estimator=:first_order, representation=:dependency_preserving),
        uncertainty_descriptions=_uncertainty_descriptions(LinearError)))
    record = (; derivative=:central, relative_step, inputs, perturbations, evaluations=count,
        observables=products)
    retained === nothing || (record=merge(record, (computations=retained,)))
    if haskey(buffers, :overlaps)
        assignment = (threshold=_MODAL_OVERLAP_THRESHOLD, frequency=axis.frequency,
            frequency_hz=frequencies(reference)[axis.frequency],
            overlap=copy(buffers.overlaps), permutation=copy(buffers.permutations))
        record = merge(record, (; mode_assignment=assignment))
    end
    return (; representation, details=ComputationDetails(record))
end

function _linear_error(method::CentralDifferences, formulation::LinearError, problem)
    Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing &&
        throw(ArgumentError("LinearError requires the Measurements extension to construct its result; " *
            "load it with `using Measurements` before compute"))
    length(problem.space) > 0 || throw(ArgumentError(
        "higher-order problem space must contain at least one core problem"))
    source_id = Commons.gridpoint_id().source_id
    populations = map(enumerate(vec(collect(points(problem.space))))) do (index, point)
        aggregate = _central_difference(method, point, formulation, problem.options, index)
        representation = Engine.retain_gridpoint(aggregate.representation,
            Commons.gridpoint_id(; source_id, problem_index=index))
        (; representation, details=aggregate.details)
    end
    return LinearErrorResult(formulation, getproperty.(populations, :representation),
        ComputationDetails(; points=getproperty.(populations, :details)))
end
