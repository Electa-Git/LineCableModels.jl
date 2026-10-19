# LineParameters computation remains independent from CableConstants.

@inline _stash!(::Nothing, ::Symbol, ::Int, ::AbstractMatrix) = nothing

@inline function _stash!(trace::NamedTuple, name::Symbol, frequency::Int, source::AbstractMatrix)
    destination = getproperty(trace, name)
    @views copyto!(destination[:, :, frequency], source)
    return nothing
end

function _solve!(
        workspace::LineParametersWorkspace{T},
        formulation::LineParametersFormulation,
        calculations::Tuple = workspace.plan.earth.calculations,
        materials::Tuple = workspace.buffers.earth.calculations
) where {T <: Real}
    input = workspace.input
    plan = workspace.plan
    buffers = workspace.buffers
    workspace.trace===nothing || empty!(workspace.trace.integrals)
    Zprimitive = buffers.Zprimitive
    Pprimitive = buffers.Pprimitive
    Zout = buffers.Zout
    Yout = buffers.Yout

    materials!(workspace, formulation)
    @debug "Starting line parameters computation"
    for frequency in 1:input.n_frequencies
        materials!(workspace, formulation, frequency, calculations, materials)
        cable_impedance!(Zprimitive, input.cable, buffers.rho_cond,
            formulation.formulas, input.jω[frequency]; workspace)
        cable_potential!(Pprimitive, input.cable, plan.shunt, buffers.dielectric_admittivity,
            input.jω[frequency], buffers.layer_coefficients, buffers.coefficients, buffers.tails)
        _stash!(workspace.trace, :Zin, frequency, Zprimitive)
        _stash!(workspace.trace, :Pin, frequency, Pprimitive)
        earth!(workspace, frequency, calculations, materials)
        _stash!(workspace.trace, :Zg, frequency, buffers.Zearth)
        _stash!(workspace.trace, :Pg, frequency, buffers.Pearth)
        impedance!(Zprimitive, workspace, frequency)
        admittance!(Pprimitive, workspace, frequency)
        reduce_line_matrices!(view(Zout, :, :, frequency), view(Yout, :, :, frequency),
            Zprimitive, Pprimitive, input.jω[frequency], plan.reduction, buffers.reduction)
    end

    return workspace
end

function _retained_details(workspace::LineParametersWorkspace{
        <:Real, <:NamedTuple, <:NamedTuple,
        <:NamedTuple, Nothing})
    # Local model diagnostics vary with geometry and selection, not the result type.
    ComputationDetails(NamedTuple{(:shunt_model,), Tuple{NamedTuple}}((workspace.plan.shunt.details,)))
end

function _retained_details(workspace::LineParametersWorkspace)
    trace = workspace.trace
    shunt = NamedTuple{(:shunt_model,), Tuple{NamedTuple}}((workspace.plan.shunt.details,))
    trace === nothing && return ComputationDetails(shunt)
    input = workspace.input
    return ComputationDetails(merge(shunt,
        (
            trace = (
            phase_map = copy(input.phase_map),
            cable_map = copy(input.cable_map),
            Zin = copy(trace.Zin),
            Pin = copy(trace.Pin),
            Zg = copy(trace.Zg),
            Pg = copy(trace.Pg),
            Z = copy(trace.Z),
            P = copy(trace.P),
            integrals = copy(trace.integrals)
        ),
        )))
end

function _finish(
        workspace::LineParametersWorkspace,
        problem::LineParametersProblem,
        formulation::LineParametersFormulation,
        basis::Units.LengthBasis;
        physical_inputs=completed_inputs(problem),
        gridpoint=Commons.gridpoint_id()
)
    impedance = copy(workspace.buffers.Zout)
    admittance = copy(workspace.buffers.Yout)
    retained = _retained_details(workspace)
    names=["cable:$(cable):$(terminal)"
           for (; cable, terminal) in problem.system.terminal_order]
    coordinates=map(workspace.plan.reduction.indices) do index
        phase=problem.system.connection_order[index]
        members=findall(==(phase), problem.system.connection_order)
        formulation.options.data.reduce_bundle && phase > 0 && length(members) > 1 ?
        "bundle:[" * join(names[members], ",") * "]" : names[index]
    end
    result = LineParameters(PhaseDomain,
        SeriesImpedance{eltype(impedance), Units.PerUnitLength}(impedance),
        ShuntAdmittance{eltype(admittance), Units.PerUnitLength}(admittance),
        workspace.input.freq,
        completion_details(merge(retained.data,
            completed_formulation(formulation, workspace),
            NamedTuple{(:inputs,:gridpoint,:coordinates),Tuple{NamedTuple,NamedTuple,Vector{String}}}(
                (physical_inputs,gridpoint,coordinates)))))
    completed = LineParameters(result, basis, workspace.input.line_length)
    all(isfinite, completed.Z.values) || throw(DomainError(completed.Z.values,
        "completed series impedance must contain only finite entries"))
    all(isfinite, completed.Y.values) || throw(DomainError(completed.Y.values,
        "completed shunt admittance must contain only finite entries"))
    return completed
end

"""
$(TYPEDSIGNATURES)

Resolve the parts of a line-parameter problem that its coaxial formulations share: the
blueprints of the designs, the line input and the completed inputs, once per problem, and
the internal potential coefficients of each formulation, whose blocks are shared where the
shunt state agrees. Return one `Resolved` per formulation.
"""
function resolve(problem::LineParametersProblem,
        formulations::Union{Tuple{LineParametersFormulation, Vararg{LineParametersFormulation}},
            AbstractVector{<:LineParametersFormulation}};
        options)
    maximum(problem.frequencies) > oftype(first(problem.frequencies), 1e8) &&
        @warn("Frequencies above 100 MHz exceed the quasi-TEM validity range.",
            max_frequency=maximum(problem.frequencies),)
    T = eltype(problem)
    blueprints = CableBlueprint{T}[DataModel.flatten(design, T)
                                   for design in problem.system.designs]
    line = lineinput(problem, blueprints)
    physical_inputs = completed_inputs(problem)
    shunts = internal_potential_coefficients(problem.system.designs, blueprints, formulations)
    return map(shunt -> Resolved(problem, (; line, shunt, physical_inputs)), shunts)
end

"""
$(TYPEDSIGNATURES)

Compute the line parameters of one point with the coaxial backend: the workspace, the
frequency loop and the completed result, from what `resolve` shares.
"""
function compute(resolved::Resolved{<:LineParametersProblem},
        formulation::LineParametersFormulation;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    execution = computation_options(LineCableModelsCoaxial, ComputationOptions(options))
    (; problem, input) = resolved
    workspace = LineParametersWorkspace(problem, formulation, execution, input.line, input.shunt)
    _solve!(workspace, formulation)
    return _finish(workspace, problem, formulation, execution.data.output_basis;
        physical_inputs = input.physical_inputs)
end

description(::LineParametersProblem; compact::Bool = false) = "Line parameters"

"""
$(TYPEDSIGNATURES)

Compute line parameters with the default formulation. See the method with a formulation.
"""
Base.@constprop :aggressive function compute(
        problem::LineParametersProblem;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions(),
        modal=nothing, modal_options::Union{NamedTuple, ComputationOptions}=ComputationOptions()
)
    return compute(problem, Formulation(); options, modal, modal_options)
end

"""
$(TYPEDSIGNATURES)

Compute the frequency-dependent line parameters of a problem or a space of problems with a
formulation or a space of formulations, through [`traverse`](@ref LineCableModels.Commons.traverse).
The formulation selects the coaxial, FEM or PSCAD backend.

The coaxial backend uses the equivalent concentric representation of the completed data
model for series impedance and ordinary radial dielectric intervals. Eligible open wire and
tape domains retain their physical geometry for the explicitly selected
`shunt_model=:boundary` calculation. The default uses annular geometry. The physical system
is normalized once into a backend-owned workspace, and all reusable numerical storage is
allocated before the frequency loop. `trace=true` copies `Zin`, `Pin`, `Zg`, `Pg`, `Z`, `P`,
`phase_map`, `cable_map` and the integration records into `details(result).data.trace`.
These arrays remain available after the workspace is reused. The result type remains
unchanged.

# Arguments

- `problem`: completed line-parameter problem, or a space of them.
- `formulation`: selected line-parameter physical methods, or a space of them.

# Keywords

- `options`: named tuple containing `verbosity`, `output_basis`, `trace`, `timing` and
  `on_result`. The optional callable `on_result(problem, index, result)` runs
  synchronously after each completed point and before the next computation. `index` is
  the point's index in the space (`1` for a scalar call). Its return value is ignored.
  Exceptions propagate. The callback must not mutate the problem or result. The default is
  `nothing`.
- `modal`: a modal-analysis formulation, a space of them or the identifier of one. The
  result is then the modal analysis of each line-parameter result.
- `modal_options`: the options of the modal analysis.

# Returns

- One [`LineParameters`](@ref) result, or a `ParametricResult` of them for a space.
"""
Base.@constprop :aggressive function compute(
        problem::ScalarOrGrid{LineParametersProblem},
        formulation::ScalarOrGrid{AbstractFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions(),
        modal=nothing, modal_options::Union{NamedTuple, ComputationOptions}=ComputationOptions()
)
    modal === nothing ||
        return compute(problem, formulation, modal; options, modal_options)
    isempty(ComputationOptions(modal_options).data) ||
        throw(ArgumentError("modal_options require a modal formulation"))
    return compute(ParametricProblem(problem, ComputationOptions(options)), formulation)
end

function computation_details(
        ::Type{<:LineParametersFormulation},
        result::LineParameters
)::ComputationDetails
    return details(result)
end

function computation_details(
        ::Type{<:LineCableModelsFEM},
        result::LineParameters
)::ComputationDetails
    return details(result)
end

function materials!(
        destination::NamedTuple,
        relation,
        model::EarthModel{T},
        frequencies::AbstractVector{T}; workspace = nothing
) where {T <: Real}
    rho, eps_r, mu_r = destination.rho, destination.eps_r, destination.mu_r
    @inbounds for row in eachindex(model.layers)
        static = EarthMaterial(model.layers[row])
        for column in eachindex(frequencies)
            material = row == firstindex(model.layers) ?
                       static :
                       constitutive(relation, static, frequencies[column]; workspace)
            rho[row, column] = material.rho
            eps_r[row, column] = material.eps_r
            mu_r[row, column] = material.mu_r
        end
    end
    return destination
end

function materials!(workspace::LineParametersWorkspace, formulation::LineParametersFormulation)
    input, buffers = workspace.input, workspace.buffers
    for (index, material) in pairs(input.cable.conductor_materials)
        buffers.rho_cond[index] = constitutive(formulation.formulas.temperature_dependence,
            material, input.temperature; workspace)
    end
    for (layer, material) in pairs(input.cable.dielectric_materials)
        buffers.rho_layer[layer] = constitutive(formulation.formulas.temperature_dependence,
            material, input.temperature; workspace)
    end
    buffers.earth.evaluated === nothing || materials!(buffers.earth.evaluated,
        formulation.formulas.earth_properties, input.earth, input.freq; workspace)
    return workspace
end

function materials!(
        workspace::LineParametersWorkspace, formulation::LineParametersFormulation,
        frequency::Int,
        calculations::Tuple = workspace.plan.earth.calculations,
        materials::Tuple = workspace.buffers.earth.calculations)
    homogenize!(workspace, frequency, formulation, calculations, materials)
    dielectric!(workspace.buffers.dielectric_admittivity, workspace.input.cable,
        workspace.buffers.rho_layer, formulation.formulas, workspace.input.freq[frequency],
        workspace.input.temperature; workspace)
    return workspace
end

@inline function _media!(destination, column::Int, air, earth)
    epsilon0 = vacuum_permittivity(typeof(earth.rho))
    mu0 = vacuum_permeability(typeof(earth.rho))
    destination.rho[1, column] = air.rho
    destination.rho[2, column] = earth.rho
    destination.epsilon[1, column] = epsilon0 * air.eps_r
    destination.epsilon[2, column] = epsilon0 * earth.eps_r
    destination.mu[1, column] = mu0 * air.mu_r
    destination.mu[2, column] = mu0 * earth.mu_r
    return nothing
end

function homogenize!(
        workspace::LineParametersWorkspace,
        frequency::Int,
        formulation::LineParametersFormulation,
        calculations::Tuple = workspace.plan.earth.calculations,
        materials::Tuple = workspace.buffers.earth.calculations
)
    foreach(calculations, materials) do calculation, destination
        homogenize!(destination, calculation, workspace, frequency,
            formulation.formulas.earth_properties)
    end
    return materials
end

# The plan decided once which earth the formula's expressions see: the layered earth or
# a reduction. The loop dispatches on that decision.
function homogenize!(destination,
        calculation::NamedTuple,
        workspace::LineParametersWorkspace, frequency_index::Int, relation)
    earth = workspace.buffers.earth
    model = workspace.input.earth
    frequency = workspace.input.freq[frequency_index]
    data = calculation.earth isa EquivalentHomogeneous.BeforeFD ? earth.static : earth.evaluated
    homogenize!(destination, calculation.earth, relation, data, model,
        calculation, frequency, frequency_index; workspace)
    return destination
end

# The layered earth reuses the layerwise FrequencyDependent values already evaluated
# during input construction. Its thicknesses were set when the buffers were allocated.
function homogenize!(
        destination,
        ::EarthModel,
        relation,
        evaluated,
        model::EarthModel,
        calculation,
        frequency,
        frequency_index::Int; workspace = nothing
)
    epsilon0=vacuum_permittivity(eltype(destination.rho))
    mu0=vacuum_permeability(eltype(destination.rho))
    # The pairs share one value object per layer, so reuse finds them identical (`===`).
    for row in axes(destination.rho, 1)
        rho=evaluated.rho[row, frequency_index]
        epsilon=epsilon0*evaluated.eps_r[row, frequency_index]
        mu=mu0*evaluated.mu_r[row, frequency_index]
        for column in eachindex(calculation.pairs)
            destination.rho[row, column]=rho
            destination.epsilon[row, column]=epsilon
            destination.mu[row, column]=mu
        end
    end
    # Every pair sees the same layered media, and the first column stands for all of them.
    isempty(calculation.pairs) || validate(@view(destination.rho[:, 1]),
        something(calculation.impedance, calculation.admittance).formula,
        @view(destination.epsilon[:, 1]), @view(destination.mu[:, 1]), destination.thickness)
    return destination
end

function homogenize!(
        destination,
        sequence::EquivalentHomogeneous.AfterFD,
        relation,
        evaluated,
        model::EarthModel,
        calculation,
        frequency,
        frequency_index::Int; workspace = nothing
)
    rho = @view evaluated.rho[:, frequency_index]
    eps_r = @view evaluated.eps_r[:, frequency_index]
    mu_r = @view evaluated.mu_r[:, frequency_index]
    air = EarthMaterial(rho[1], eps_r[1], mu_r[1])
    foreach(calculation.reductions) do reduction
        @inbounds for column in reduction.pairs
            functor = Functor(sequence.rule, (; rho, eps_r, mu_r, model,
                pair = calculation.pairs[column].physical, frequency, reduction.options);
                workspace)
            earth = validate(reduction.expression(functor, workspace), sequence.rule)
            _media!(destination, column, air, earth)
            validate(@view(destination.rho[:, column]),
                something(calculation.impedance, calculation.admittance).formula,
                @view(destination.epsilon[:, column]), @view(destination.mu[:, column]),
                destination.thickness)
        end
    end
    return destination
end

function homogenize!(
        destination,
        sequence::EquivalentHomogeneous.BeforeFD,
        relation,
        static,
        model::EarthModel,
        calculation,
        frequency,
        frequency_index::Int; workspace = nothing
)
    air = EarthMaterial(static.rho[1], static.eps_r[1], static.mu_r[1])
    foreach(calculation.reductions) do reduction
        @inbounds for column in reduction.pairs
            functor = Functor(sequence.rule, (; static.rho, static.eps_r, static.mu_r, model,
                pair = calculation.pairs[column].physical, frequency, reduction.options);
                workspace)
            reconstructed = validate(reduction.expression(functor, workspace), sequence.rule)
            earth = constitutive(relation, reconstructed, frequency; workspace)
            _media!(destination, column, air, earth)
            validate(@view(destination.rho[:, column]),
                something(calculation.impedance, calculation.admittance).formula,
                @view(destination.epsilon[:, column]), @view(destination.mu[:, column]),
                destination.thickness)
        end
    end
    return destination
end
