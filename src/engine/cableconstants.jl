"""
$(TYPEDEF)

Store earth-free cable constants per unit length for one or more independent
concentric assemblies.

Every entry of `cores`, `R`, `L`, `C`, and `G` describes one assembly. Values
are evaluated at `frequency`; `R`, `L`, `C`, and `G` use Ω/m, H/m, F/m, and
S/m respectively.

$(TYPEDFIELDS)
"""
struct CableConstants{T <: Real, D <: ComputationDetails} <: AbstractCoreResult
    "Innermost active terminal of each concentric assembly."
    cores::Vector{Symbol}
    "Series resistance per unit length [Ω/m]."
    R::Vector{T}
    "Series inductance per unit length [H/m]."
    L::Vector{T}
    "Shunt capacitance per unit length [F/m]."
    C::Vector{T}
    "Shunt conductance per unit length [S/m]."
    G::Vector{T}
    "Evaluation frequency [Hz]."
    frequency::T

    "Requested/resolved formulas and blueprint construction diagnostics."
    details::D

    function CableConstants{T}(
            cores::Vector{Symbol},
            R::Vector{T},
            L::Vector{T},
            C::Vector{T},
            G::Vector{T},
            frequency::T,
            details::ComputationDetails = ComputationDetails()
    ) where {T <: Real}
        count = length(cores)
        iszero(count) && throw(ArgumentError(
            "cable constants require at least one concentric assembly",
        ))
        all(length(values) == count for values in (R, L, C, G)) ||
            throw(DimensionMismatch(
                "cores, R, L, C, and G must contain the same number of entries",
            ))
        allunique(cores) || throw(ArgumentError(
            "cable-constant core names must be unique",
        ))
        isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
            frequency,
            "cable-constant frequency must be positive and finite"
        ))
        all(isfinite, Iterators.flatten((R, L, C, G))) || throw(DomainError(
            (R, L, C, G),
            "cable constants must be finite"
        ))
        return new{T, typeof(details)}(cores, R, L, C, G, frequency, details)
    end
end

function CableConstants(
        cores::AbstractVector{Symbol},
        R::AbstractVector{<:Real},
        L::AbstractVector{<:Real},
        C::AbstractVector{<:Real},
        G::AbstractVector{<:Real},
        frequency::Real,
        details::ComputationDetails = ComputationDetails()
)
    T = promote_type(
        eltype(R), eltype(L), eltype(C), eltype(G), typeof(float(frequency))
    )
    return CableConstants{T}(
        collect(Symbol, cores),
        T.(R),
        T.(L),
        T.(C),
        T.(G),
        convert(T, float(frequency)),
        details
    )
end

function CableConstants(
        R::Real,
        L::Real,
        C::Real,
        G::Real = 0;
        core::Symbol = :core,
        frequency::Real = 50
)
    values = promote(R, L, C, G, frequency)
    T = typeof(first(values))
    return CableConstants{T}(
        Symbol[core],
        T[values[1]],
        T[values[2]],
        T[values[3]],
        T[values[4]],
        values[5]
    )
end

function Base.:(==)(left::CableConstants, right::CableConstants)
    left.cores == right.cores && left.R == right.R && left.L == right.L &&
        left.C == right.C && left.G == right.G && left.frequency == right.frequency
end

Base.length(constants::CableConstants) = length(constants.cores)
details(constants::CableConstants) = constants.details
Base.size(constants::CableConstants) = (length(constants),)
function Base.eltype(::Type{CableConstants{T}}) where {T}
    NamedTuple{
        (:core, :R, :L, :C, :G),
        Tuple{Symbol, T, T, T, T}
    }
end
Base.firstindex(constants::CableConstants) = firstindex(constants.cores)
Base.lastindex(constants::CableConstants) = lastindex(constants.cores)

function Base.getindex(constants::CableConstants, index::Integer)
    return (
        core = constants.cores[index],
        R = constants.R[index],
        L = constants.L[index],
        C = constants.C[index],
        G = constants.G[index]
    )
end

function Base.iterate(constants::CableConstants, state::Int = 1)
    state > length(constants) && return nothing
    return constants[state], state + 1
end

observe(constants::CableConstants, ::typeof(R)) = constants.R
observe(constants::CableConstants, ::typeof(L)) = constants.L
observe(constants::CableConstants, ::typeof(C)) = constants.C
observe(constants::CableConstants, ::typeof(G)) = constants.G

R(constants::CableConstants) = observe(constants, R)
L(constants::CableConstants) = observe(constants, L)
C(constants::CableConstants) = observe(constants, C)
G(constants::CableConstants) = observe(constants, G)
basis(::CableConstants) = :pul
resistance(constants::CableConstants) = observe(constants, R)
inductance(constants::CableConstants) = observe(constants, L)
capacitance(constants::CableConstants) = observe(constants, C)
conductance(constants::CableConstants) = observe(constants, G)
observables(::Type{<:CableConstants}) = (R, L, C, G)

function publication_table(
        source::CableConstants,
        requests::Tuple,
        observations::Tuple,
        ::NamedTuple
)
    names = map(payload -> Symbol(Units.symbol(payload.quantity)), observations)
    length(unique(names)) == length(names) || throw(ArgumentError(
        "cable-constant publication quantities must be distinct",
    ))
    all(payload -> length(payload.values) == length(source), observations) ||
        throw(DimensionMismatch(
            "cable-constant observations must align with the assembly count",
        ))
    values = NamedTuple{names}(map(payload -> collect(payload.values), observations))
    contract = NamedTuple{names}(map(observations) do payload
        (; quantity = payload.quantity, unit = payload.unit)
    end)
    return (
        columns = merge((core = copy(source.cores),), values),
        row_order = (:core, names...),
        observation_columns = contract
    )
end

"""
$(TYPEDEF)

Define one earth-free cable-constant calculation for a completed cable design.

The innermost terminal of every concentric assembly is active and every
additional outward terminal, when present, is grounded. The calculation is
restricted to the 50 Hz or 60 Hz base frequency used by cable datasheets.

$(TYPEDFIELDS)
"""
struct CableConstantsProblem{
    T <: Real,
    D <: CableDesign
} <: AbstractProblemDefinition
    "Completed physical cable design."
    design::D
    "Operating temperature [°C]."
    temperature::T
    "Evaluation frequency [Hz]."
    frequency::T

    function CableConstantsProblem{T, D}(
            design::D,
            temperature::T,
            frequency::T
    ) where {T <: Real, D <: CableDesign}
        return validate(new{T, D}(design, temperature, frequency))
    end
end

Base.eltype(::CableConstantsProblem{T}) where {T} = T
Base.eltype(::Type{<:CableConstantsProblem{T}}) where {T} = T

function validate(problem::CableConstantsProblem)
    validate(problem.design)
    isfinite(problem.temperature) || throw(DomainError(
        problem.temperature,
        "cable-constant temperature must be finite"
    ))
    problem.frequency in (oftype(problem.frequency, 50), oftype(problem.frequency, 60)) ||
        throw(DomainError(
            problem.frequency,
            "cable-constant base frequency must be 50 Hz or 60 Hz"
        ))
    return problem
end

"""
$(TYPEDSIGNATURES)

Construct an earth-free cable-constant problem.

# Keywords

- `temperature=20`: Operating temperature [°C].
- `frequency=50`: Base frequency, either 50 Hz or 60 Hz.

# Returns

- A validated [`CableConstantsProblem`](@ref).
"""
function CableConstantsProblem(
        design::CableDesign;
        temperature::Real = 20,
        frequency::Real = 50
)
    T = promote_type(
        eltype(design), typeof(float(temperature)), typeof(float(frequency))
    )
    value = convert(T, float(frequency))
    return CableConstantsProblem{T, typeof(design)}(
        design,
        convert(T, float(temperature)),
        value
    )
end

"""
$(TYPEDEF)

Select the conductor and dielectric formulas used by a cable-constant
calculation.

$(TYPEDFIELDS)
"""
struct CableConstantsFormulation{
    M <: NamedTuple,
    O <: FormulationOptions,
    D <: NamedTuple
} <: AbstractFormulation
    "Registered physical formula selections."
    methods::M
    "Cable-constant formulation options."
    options::O
    "Requested declarations retained before formula-owner resolution."
    definitions::D
end

function formulation_options(
        ::Type{CableConstantsFormulation},
        options::FormulationOptions
)::FormulationOptions
    isempty(options.data) || throw(ArgumentError(
        "unknown cable-constant formulation options: $(keys(options.data))"))
    return FormulationOptions()
end

description(::Type{<:CableConstantsFormulation}; compact::Bool = false) = "Cable constants"
function description(::CableConstantsFormulation; compact::Bool = false)
    description(CableConstantsFormulation; compact)
end
formula_id(::Type{<:CableConstantsFormulation}) = :cable_constants
formula_id(::CableConstantsFormulation) = :cable_constants
formulation_options(value::CableConstantsFormulation) = value.options
function Base.pairs(::Type{CableConstantsFormulation}; quantity = nothing)
    return pairs((;
        (key=>family
    for (key, family) in pairs(LineParametersFormulation; quantity)
    if key ∉ (:earth_impedance, :earth_admittance, :earth_properties))...))
end
function description(::Type{CableConstantsFormulation}, slot::Val)
    description(LineParametersFormulation, slot)
end
function Base.pairs(value::CableConstantsFormulation; quantity = nothing)
    pairs(CableConstantsFormulation,
        (methods = value.methods, requested = value.definitions,
            options = value.options.data); quantity)
end
function Base.pairs(::Type{CableConstantsFormulation}, retained::NamedTuple; quantity = nothing)
    pairs(LineParametersFormulation, retained; quantity, owner = CableConstantsFormulation)
end

"""Expose requested and resolved local formulas for passive scientific records."""
function Base.NamedTuple(value::CableConstantsFormulation)
    record = function (selected)
        selected === nothing && return nothing
        selected isa Symbol && return NamedTuple(formula(selected))
        selected isa NamedTuple && return map(record, selected)
        return NamedTuple(selected)
    end
    Record = NamedTuple{(:backend, :requested, :methods, :options),
        Tuple{Symbol, NamedTuple, NamedTuple, NamedTuple}}
    return Record((:cable_constants, map(record, value.definitions),
        map(record, value.methods), value.options.data))
end

function _constants_formulation(
        internal_impedance,
        insulation_impedance,
        shunt_model,
        insulation_admittance,
        semicon_admittance,
        pipe_impedance,
        temperature_dependence,
        options::FormulationOptions
)
    methods = (
        internal_impedance = Formulation(InternalImpedance.Formula, internal_impedance),
        insulation_impedance = InsulationImpedance.Formula(insulation_impedance),
        shunt_model = ShuntModel.Formula(shunt_model),
        insulation_admittance = InsulationAdmittance.Formula(insulation_admittance),
        semicon_admittance = SemiconAdmittance.Formula(semicon_admittance),
        pipe_impedance = PipeImpedance.Formula(pipe_impedance),
        temperature_dependence = temperature_dependence === nothing ? nothing :
                                 TemperatureDependent.Formula(temperature_dependence)
    )
    internal_impedance isa NamedTuple &&
        (internal_impedance = NamedTuple{keys(methods.internal_impedance)}(internal_impedance))
    return CableConstantsFormulation(
        methods,
        formulation_options(CableConstantsFormulation, options),
        (; internal_impedance, insulation_impedance, shunt_model, insulation_admittance,
            semicon_admittance, pipe_impedance, temperature_dependence)
    )
end

"""
$(TYPEDSIGNATURES)

Construct the cable-constant formula bundle.

Each formula slot and the complete `options` tuple accepts a scalar selection
or an explicit [`Grid`](@ref LineCableModels.ParametricBuilder.Grid)/
[`Gridspace`](@ref LineCableModels.ParametricBuilder.Gridspace) source. Scalar
inputs return one [`CableConstantsFormulation`](@ref); varying inputs return a
`Gridspace{CableConstantsFormulation}` of completed formulations.

# Keywords

- `internal_impedance`: Conductor surface-impedance recipe.
- `insulation_impedance`: Longitudinal insulation-impedance recipe.
- `shunt_model`: Local geometry model; `:default`/`:coaxial` uses annuli,
  `:boundary` explicitly prepares lossless open-screen coupling.
- `insulation_admittance`: Insulation constitutive relation.
- `semicon_admittance`: Semiconducting-layer constitutive relation.
- `pipe_impedance`: Pipe-type selection; the coaxial pipe implementation is not
  yet available. Ordinary concentric assemblies have no additional pipe term.
- `temperature_dependence`: Resistivity law; `:default` is linear and `nothing`
  retains reference resistivity.
- `options`: Formulation controls; currently empty.
- `combine`: `:product` or `:zip` composition among varying fields.
"""
function CableConstantsFormulation(;
        internal_impedance = formula(:default),
        insulation_impedance = formula(:default),
        shunt_model = formula(:default),
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
        pipe_impedance = formula(:default),
        temperature_dependence = formula(:default),
        options = FormulationOptions(),
        combine::Symbol = :product
)
    values = (
        internal_impedance,
        insulation_impedance,
        shunt_model,
        insulation_admittance,
        semicon_admittance,
        pipe_impedance,
        temperature_dependence,
        options
    )
    return parameterize(
        CableConstantsFormulation,
        (inputs...) -> _constants_formulation(inputs[1:(end - 1)]...,
            last(inputs) isa NamedTuple ? FormulationOptions(last(inputs)) : last(inputs)),
        values;
        combine
    )
end

function computation_options(
        ::Type{<:CableConstantsFormulation},
        options::ComputationOptions
)::ComputationOptions
    isempty(options.data) || throw(ArgumentError(
        "CableConstants compute does not accept computation options",
    ))
    return ComputationOptions()
end

"""
$(TYPEDEF)

Own the local cable arrays, corrected resistivities, and reusable matrices for
one cable-constant calculation. The constructor consumes the completed
blueprint without retaining a duplicate representation.

$(TYPEDFIELDS)
"""
struct CableConstantsWorkspace{T <: Real, L, B}
    "Concrete array payload used by local primitive assemblers."
    cable::L
    "Temperature-corrected conductor resistivities [Ω·m]."
    rho::Vector{T}
    "Reusable primitive, reduction, and result storage."
    buffers::B
end

function CableConstantsWorkspace(
        problem::CableConstantsProblem{T},
        formulation::CableConstantsFormulation,
        blueprint::CableBlueprint{T}
) where {T <: Real}
    return CableConstantsWorkspace(
        problem,
        formulation,
        LocalCableData(blueprint)
    )
end

function CableConstantsWorkspace(
        problem::CableConstantsProblem{T},
        formulation::CableConstantsFormulation,
        cable::LocalCableData{T}
) where {T <: Real}
    @inbounds for assembly in cable.assemblies
        isempty(cable.dielectric_ranges[first(assembly)]) && throw(ArgumentError(
            "assembly core :$(cable.terminals[first(assembly)]) has no radial dielectric path",
        ))
    end
    count = length(cable.terminals)
    rho = Vector{T}(undef, length(cable.conductor_materials))

    maximum_size = maximum(length, cable.assemblies)
    removed = maximum_size - 1
    buffers = (
        Z = Matrix{Complex{T}}(undef, count, count),
        Y = Matrix{Complex{T}}(undef, count, count),
        reduced = Matrix{Complex{T}}(undef, 1, 1),
        factor = Matrix{Complex{T}}(undef, removed, removed),
        coupling = Matrix{Complex{T}}(undef, 1, removed),
        right_hand_side = Matrix{Complex{T}}(undef, removed, 1),
        indices = collect(1:maximum_size),
        layer_coefficients = Vector{Complex{T}}(
            undef, length(cable.dielectric_materials)
        ),
        dielectric_admittivity = Vector{Complex{T}}(undef, length(cable.dielectric_materials)),
        R = Vector{T}(undef, length(cable.assemblies)),
        L = Vector{T}(undef, length(cable.assemblies)),
        C = Vector{T}(undef, length(cable.assemblies)),
        G = Vector{T}(undef, length(cable.assemblies)),
        quadrature = integration_workspace(typeof(float(nominal(one(T)))), Complex{T}; size = 0),
        observations = nothing
    )
    selected_internal = formulation.methods.internal_impedance
    internal = if selected_internal isa NamedTuple
        kinds = any(>(0), cable.r_in) ? (:inner, :outer, :transfer) : (:outer,)
        Tuple(unique(Formulation(selected_internal, Val(kind)) for kind in kinds))
    else
        selected_internal
    end
    buffers = initialize_buffers(merge(formulation.methods, (internal_impedance = internal,)),
        T, cable, (;), buffers)
    return CableConstantsWorkspace{T, typeof(cable), typeof(buffers)}(
        cable, rho, buffers
    )
end

function _solve!(
        workspace::CableConstantsWorkspace{T},
        problem::CableConstantsProblem{T},
        formulation::CableConstantsFormulation
) where {T <: Real}
    buffers = workspace.buffers
    ω = 2 * (one(T) * π) * problem.frequency
    s = complex(zero(T), ω)
    for (index, material) in pairs(workspace.cable.conductor_materials)
        workspace.rho[index] = constitutive(formulation.methods.temperature_dependence,
            material, problem.temperature; workspace)
    end
    dielectric!(buffers.dielectric_admittivity, workspace.cable, formulation.methods,
        problem.frequency, problem.temperature; workspace)
    cable_impedance!(
        buffers.Z,
        workspace.cable,
        workspace.rho,
        formulation.methods,
        s; workspace
    )
    cable_admittance!(
        buffers.Y,
        workspace.cable,
        buffers.dielectric_admittivity,
        s,
        buffers.layer_coefficients
    )
    keep = @view buffers.indices[1:1]
    @inbounds for (assembly, chain) in pairs(workspace.cable.assemblies)
        count = length(chain)
        Z = @view buffers.Z[chain, chain]
        eliminate = @view buffers.indices[2:count]
        factor = @view buffers.factor[1:(count - 1), 1:(count - 1)]
        coupling = @view buffers.coupling[:, 1:(count - 1)]
        right_hand_side = @view buffers.right_hand_side[1:(count - 1), :]
        kronify!(
            Z,
            keep,
            eliminate,
            buffers.reduced,
            factor,
            coupling,
            right_hand_side
        )
        equivalent = buffers.reduced[1, 1]
        buffers.R[assembly] = real(equivalent)
        buffers.L[assembly] = imag(equivalent) / ω

        Y = buffers.Y[first(chain), first(chain)]
        iszero(Y) && throw(ArgumentError(
            "assembly core :$(workspace.cable.terminals[first(chain)]) has no finite radial dielectric path",
        ))
        buffers.G[assembly] = real(Y)
        buffers.C[assembly] = imag(Y) / ω
    end
    return CableConstants(
        Symbol[workspace.cable.terminals[first(chain)]
               for chain in workspace.cable.assemblies],
        buffers.R,
        buffers.L,
        buffers.C,
        buffers.G,
        problem.frequency,
        ComputationDetails(NamedTuple{(:shunt_model, :formulations),
            Tuple{NamedTuple, NamedTuple}}((workspace.cable.shunt_details, NamedTuple(formulation))))
    )
end

"""
$(TYPEDSIGNATURES)

Compute earth-free cable constants with the default coaxial formulation.
"""
function compute(
        problem::CableConstantsProblem;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(
        LineCableModelsCoaxial(),
        problem,
        CableConstantsFormulation();
        options
    )
end

"""
$(TYPEDSIGNATURES)

Compute earth-free cable constants with an explicit formula bundle.
"""
function compute(
        problem::CableConstantsProblem,
        formulation::CableConstantsFormulation;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(LineCableModelsCoaxial(), problem, formulation; options)
end

function compute(
        problem::CableConstantsProblem,
        formulations::AbstractVector{<:CableConstantsFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return compute(LineCableModelsCoaxial(), problem, formulations; options)
end

"""
$(TYPEDSIGNATURES)

Compute earth-free cable constants through the coaxial backend tag.
"""
function compute(
        engine::LineCableModelsCoaxial,
        problem::CableConstantsProblem,
        formulation::CableConstantsFormulation = CableConstantsFormulation();
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    values = compute(
        engine,
        problem,
        typeof(formulation)[formulation];
        options
    )
    return first(values)
end

function compute(
        engine::LineCableModelsCoaxial,
        problem::CableConstantsProblem,
        formulations::AbstractVector{<:CableConstantsFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    options = options isa NamedTuple ? ComputationOptions(options) : options
    computation_options(CableConstantsFormulation, options)
    isempty(formulations) && throw(ArgumentError(
        "cable-constant formulation collections cannot be empty",
    ))
    validate(problem)
    for formulation in formulations
        Formulation(engine, formulation.methods.pipe_impedance, problem.design)
    end
    blueprints = flatten(engine, [problem.design], eltype(problem), formulations)
    cables = [LocalCableData(first(blueprints))]
    for index in 2:length(blueprints)
        previous = findfirst(other -> other === blueprints[index], blueprints)
        push!(cables, previous < index ? cables[previous] :
                      LocalCableData(blueprints[index]))
    end
    return map(formulations, cables) do formulation, cable
        workspace = CableConstantsWorkspace(problem, formulation, cable)
        _solve!(workspace, problem, formulation)
    end
end

"""
$(TYPEDSIGNATURES)

Calculate earth-free constants directly from one completed cable design.
"""
function CableConstants(
        design::CableDesign;
        temperature::Real = 20,
        frequency::Real = 50,
        formulation::CableConstantsFormulation = CableConstantsFormulation(),
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions()
)
    problem = CableConstantsProblem(design; temperature, frequency)
    return compute(problem, formulation; options)
end

function computation_details(
        ::Type{<:CableConstantsFormulation},
        result::CableConstants
)::ComputationDetails
    return details(result)
end
