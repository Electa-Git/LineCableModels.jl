"""
$(TYPEDEF)

Store earth-free cable constants per unit length for one or more independent
concentric assemblies.

Every entry of `cores`, `R`, `L`, `C`, and `G` describes one assembly. Values
are evaluated at `frequency`; `R`, `L`, `C`, and `G` use Ω/m, H/m, F/m, and
S/m respectively.

$(TYPEDFIELDS)
"""
struct CableConstants{T <: Real} <: AbstractCoreResult
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

    function CableConstants{T}(
            cores::Vector{Symbol},
            R::Vector{T},
            L::Vector{T},
            C::Vector{T},
            G::Vector{T},
            frequency::T
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
            "cable-constant frequency must be positive and finite",
        ))
        all(isfinite, Iterators.flatten((R, L, C, G))) || throw(DomainError(
            (R, L, C, G),
            "cable constants must be finite",
        ))
        return new{T}(cores, R, L, C, G, frequency)
    end
end

function CableConstants(
        cores::AbstractVector{Symbol},
        R::AbstractVector{<:Real},
        L::AbstractVector{<:Real},
        C::AbstractVector{<:Real},
        G::AbstractVector{<:Real},
        frequency::Real
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
        convert(T, float(frequency))
    )
end

function CableConstants(
        R::Real,
        L::Real,
        C::Real,
        G::Real = 0;
        core::Symbol = :core,
        frequency::Real = 1e-3
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

Base.:(==)(left::CableConstants, right::CableConstants) =
    left.cores == right.cores && left.R == right.R && left.L == right.L &&
    left.C == right.C && left.G == right.G && left.frequency == right.frequency

Base.length(constants::CableConstants) = length(constants.cores)

"Return the only assembly as a scalar named tuple."
function Base.only(constants::CableConstants)
    length(constants) == 1 || throw(ArgumentError(
        "CableConstants contains $(length(constants)) assemblies; expected exactly one",
    ))
    return (
        core = only(constants.cores),
        R = only(constants.R),
        L = only(constants.L),
        C = only(constants.C),
        G = only(constants.G)
    )
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
resistance(constants::CableConstants) = R(constants)
inductance(constants::CableConstants) = L(constants)
capacitance(constants::CableConstants) = C(constants)
conductance(constants::CableConstants) = G(constants)
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

function _chains(components)
    isempty(components) && throw(ArgumentError(
        "cable constants require at least one flattened component",
    ))
    starts = Int[1]
    positions = [first(components).conductor.position]
    @inbounds for index in 2:length(components)
        position = components[index].conductor.position
        if !DataModel.same_radial_position(position, last(positions))
            any(reference -> DataModel.same_radial_position(position, reference), positions) &&
                throw(ArgumentError(
                    "a concentric assembly cannot reappear after another assembly",
                ))
            push!(starts, index)
            push!(positions, position)
        end
    end
    return UnitRange{Int}[
        start:(index == length(starts) ? length(components) : starts[index + 1] - 1)
        for (index, start) in pairs(starts)
    ]
end

"""
$(TYPEDEF)

Define one earth-free cable-constant calculation for a completed cable design.

The internally constructed system contains one design at the origin. The
innermost terminal of every concentric assembly is active and every outward
terminal is grounded.

$(TYPEDFIELDS)
"""
struct CableConstantsProblem{
    T <: Real,
    S <: LineCableSystem
} <: AbstractProblemDefinition
    "Neutral one-design cable system."
    system::S
    "Operating temperature [°C]."
    temperature::T
    "Evaluation frequency [Hz]."
    frequency::T
end

Base.eltype(::CableConstantsProblem{T}) where {T} = T
Base.eltype(::Type{<:CableConstantsProblem{T}}) where {T} = T

function _check_cable_constants_problem(problem::CableConstantsProblem)
    system = problem.system
    validate(system)
    length(system.designs) == 1 || throw(ArgumentError(
        "a cable-constant problem requires exactly one cable design",
    ))
    system.environment === nothing || throw(ArgumentError(
        "a cable-constant problem cannot define an external environment",
    ))
    isone(system.line_length) || throw(ArgumentError(
        "a cable-constant problem uses a fixed one-metre reference length",
    ))
    pose = only(system.positions)
    all(iszero, (pose.x, pose.y, pose.φ)) || throw(ArgumentError(
        "a cable-constant problem requires its design at the origin",
    ))
    isfinite(problem.temperature) || throw(DomainError(
        problem.temperature,
        "cable-constant temperature must be finite",
    ))
    isfinite(problem.frequency) && problem.frequency > zero(problem.frequency) ||
        throw(DomainError(
            problem.frequency,
            "cable-constant frequency must be positive and finite",
        ))

    design = only(system.designs)
    components = DataModel.flatten(design, problem.frequency, eltype(problem))
    chains = _chains(components)
    all(length(chain) >= 2 for chain in chains) || throw(ArgumentError(
        "every concentric assembly requires an outward grounded terminal",
    ))
    names = getproperty.(components, :name)
    names == design.terminal_order || throw(DimensionMismatch(
        "DataModel terminal order differs from canonical flattened components",
    ))
    system.terminal_order == [(cable = 1, terminal = name) for name in names] ||
        throw(DimensionMismatch(
            "cable-system terminal order differs from canonical flattened components",
        ))
    expected = zeros(Int, length(components))
    for (assembly, chain) in pairs(chains)
        expected[first(chain)] = assembly
        isempty(first(components[chain]).dielectric.layers) && throw(ArgumentError(
            "assembly core :$(first(components[chain]).name) has no radial dielectric path to its grounded terminal",
        ))
    end
    system.connection_order == expected || throw(DimensionMismatch(
        "cable-constant terminal connections must retain only each assembly core",
    ))
    for region in design.geometry.regions
        DataModel.temperature_factor(
            region.source.material.alpha,
            problem.temperature,
            region.source.material.T0
        )
    end
    return nothing
end

function Validation.rules(::Type{<:CableConstantsProblem})
    (Validation.OwnerRule(:cable_constants_problem, _check_cable_constants_problem),)
end

"""
$(TYPEDSIGNATURES)

Construct an earth-free cable-constant problem.

# Keywords

- `temperature=20`: Operating temperature [°C].
- `frequency=1e-3`: Positive evaluation frequency [Hz].

# Returns

- A validated [`CableConstantsProblem`](@ref).
"""
function CableConstantsProblem(
        design::CableDesign;
        temperature::Real = 20,
        frequency::Real = 1e-3
)
    T = promote_type(
        eltype(design), typeof(float(temperature)), typeof(float(frequency))
    )
    value = convert(T, float(frequency))
    components = DataModel.flatten(design, value, T)
    chains = _chains(components)
    connections = zeros(Int, length(components))
    @inbounds for (assembly, chain) in pairs(chains)
        connections[first(chain)] = assembly
    end
    system = build(
        LineCableSystem,
        design,
        DataModel.Pose2(zero(T), zero(T), zero(T)),
        connections,
        nothing,
        "$(design.cable_id)_constants",
        one(T)
    )
    problem = CableConstantsProblem{T, typeof(system)}(
        system,
        convert(T, float(temperature)),
        value
    )
    return validate(problem)
end

"""
$(TYPEDEF)

Select the conductor and dielectric formulas used by a cable-constant
calculation.

$(TYPEDFIELDS)
"""
struct CableConstantsFormulation{
    M <: NamedTuple,
    O <: NamedTuple
} <: AbstractFormulation
    "Registered physical formula selections."
    methods::M
    "Cable-constant formulation options."
    options::O
end

function formulation_options(
        ::Val{CableConstantsFormulation},
        options::NamedTuple
)::FormulationOptions
    allowed = (:temperature_correction,)
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown cable-constant formulation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge((temperature_correction = true,), options)
    normalized.temperature_correction isa Bool || throw(ArgumentError(
        "temperature_correction must be Bool",
    ))
    return (temperature_correction = normalized.temperature_correction,)
end

"""
$(TYPEDSIGNATURES)

Construct the cable-constant formula bundle.

# Keywords

- `internal_impedance`: Conductor surface-impedance recipe.
- `insulation_impedance`: Longitudinal insulation-impedance recipe.
- `insulation_admittance`: Insulation constitutive relation.
- `semicon_admittance`: Semiconducting-layer constitutive relation.
- `options`: Cable-constant formulation options.
"""
function CableConstantsFormulation(;
        internal_impedance = formula(:Schelkunoff1934),
        insulation_impedance = formula(:Ametani1980),
        insulation_admittance = formula(:Marti2001),
        semicon_admittance = formula(:Ametani2004),
        options::NamedTuple = (;)
)
    methods = (
        internal_impedance = _internal_impedance_formula(internal_impedance),
        insulation_impedance = _insulation_impedance_formula(insulation_impedance),
        insulation_admittance = _insulation_admittance_formula(insulation_admittance),
        semicon_admittance = _semicon_admittance_formula(semicon_admittance)
    )
    return CableConstantsFormulation(
        methods,
        formulation_options(Val(CableConstantsFormulation), options)
    )
end

function computation_options(
        ::Val{CableConstantsProblem},
        options::NamedTuple
)::ComputationOptions
    isempty(options) || throw(ArgumentError(
        "CableConstants compute does not accept computation options",
    ))
    return (;)
end

struct CableConstantsWorkspace{T <: Real, V, B}
    components::V
    chains::Vector{UnitRange{Int}}
    rho::Vector{T}
    mu_ins::Vector{T}
    buffers::B
end

function CableConstantsWorkspace(
        problem::CableConstantsProblem{T},
        formulation::CableConstantsFormulation
) where {T <: Real}
    design = only(problem.system.designs)
    components = DataModel.flatten(design, problem.frequency, T)
    chains = _chains(components)
    count = length(components)
    rho = Vector{T}(undef, count)
    mu_ins = Vector{T}(undef, count)
    @inbounds for index in eachindex(components)
        component = components[index]
        material = component.conductor.material
        rho[index] = material.rho
        if formulation.options.temperature_correction
            rho[index] *= DataModel.temperature_factor(
                material.alpha,
                problem.temperature,
                material.T0
            )
        end
        layers = component.dielectric.layers
        mu_ins[index] = isempty(layers) ? one(T) :
                        DataModel.equivalent_dielectric_permeability(
            layers,
            component.conductor.num_turns,
            component.conductor.r_ex,
            component.dielectric.r_ex
        )
    end

    maximum_size = maximum(length, chains)
    removed = maximum_size - 1
    buffers = (
        Z = Matrix{Complex{T}}(undef, maximum_size, maximum_size),
        reduced = Matrix{Complex{T}}(undef, 1, 1),
        factor = Matrix{Complex{T}}(undef, removed, removed),
        coupling = Matrix{Complex{T}}(undef, 1, removed),
        right_hand_side = Matrix{Complex{T}}(undef, removed, 1),
        indices = collect(1:maximum_size),
        R = Vector{T}(undef, length(chains)),
        L = Vector{T}(undef, length(chains)),
        C = Vector{T}(undef, length(chains)),
        G = Vector{T}(undef, length(chains))
    )
    return CableConstantsWorkspace{T, typeof(components), typeof(buffers)}(
        components, chains, rho, mu_ins, buffers
    )
end

function _series!(
        destination::AbstractMatrix{Complex{T}},
        components,
        chain::UnitRange{Int},
        rho::AbstractVector{T},
        mu_ins::AbstractVector{T},
        formulation::CableConstantsFormulation,
        s::Complex{T}
) where {T <: Real}
    fill!(destination, zero(Complex{T}))
    count = length(chain)
    @inbounds for position in count:-1:1
        index = chain[position]
        component = components[index]
        conductor = component.conductor
        interaction = formulation.methods.internal_impedance(
            conductor.r_in,
            conductor.r_ex,
            rho[index],
            conductor.material.mu_r,
            s
        )
        outside = interaction(Val(:outer))
        inside = if position < count
            next_index = chain[position + 1]
            next_conductor = components[next_index].conductor
            formulation.methods.internal_impedance(
                next_conductor.r_in,
                next_conductor.r_ex,
                rho[next_index],
                next_conductor.material.mu_r,
                s
            )(Val(:inner))
        else
            zero(outside)
        end
        mutual = interaction(Val(:mutual))
        insulation = formulation.methods.insulation_impedance(
            conductor.r_ex,
            component.dielectric.r_ex,
            mu_ins[index],
            s
        )
        loop = outside + inside + insulation
        if position > 1
            for row in 1:(position - 1), column in 1:(position - 1)
                destination[row, column] += loop - 2mutual
            end
            for row in 1:(position - 1)
                destination[position, row] += loop - mutual
                destination[row, position] += loop - mutual
            end
        end
        destination[position, position] += loop
    end
    return destination
end

function _shunt(
        component,
        problem::CableConstantsProblem{T},
        formulation::CableConstantsFormulation,
        s::Complex{T}
) where {T <: Real}
    radial = zero(Complex{T})
    for layer in component.dielectric.layers
        selected = if layer.material.kind === :insulator
            formulation.methods.insulation_admittance
        elseif layer.material.kind === :semicon
            formulation.methods.semicon_admittance
        else
            throw(ArgumentError(
                "cable-constant dielectric layers must be :insulator or :semicon; got :$(layer.material.kind)",
            ))
        end
        material = constitutive(
            selected,
            layer.material,
            problem.frequency,
            problem.temperature
        )
        radial += inv(layer_admittance(
            layer.r_in,
            layer.r_ex,
            admittivity(material, s)
        ))
    end
    iszero(radial) && throw(ArgumentError(
        "a cable-constant assembly requires a finite radial dielectric path",
    ))
    return inv(radial)
end

function _solve!(
        workspace::CableConstantsWorkspace{T},
        problem::CableConstantsProblem{T},
        formulation::CableConstantsFormulation
) where {T <: Real}
    buffers = workspace.buffers
    ω = 2 * (one(T) * π) * problem.frequency
    s = complex(zero(T), ω)
    keep = @view buffers.indices[1:1]
    @inbounds for (assembly, chain) in pairs(workspace.chains)
        count = length(chain)
        Z = @view buffers.Z[1:count, 1:count]
        _series!(
            Z,
            workspace.components,
            chain,
            workspace.rho,
            workspace.mu_ins,
            formulation,
            s
        )
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

        Y = _shunt(
            workspace.components[first(chain)],
            problem,
            formulation,
            s
        )
        buffers.G[assembly] = real(Y)
        buffers.C[assembly] = imag(Y) / ω
    end
    return CableConstants(
        Symbol[workspace.components[first(chain)].name for chain in workspace.chains],
        buffers.R,
        buffers.L,
        buffers.C,
        buffers.G,
        problem.frequency
    )
end

"""
$(TYPEDSIGNATURES)

Compute earth-free cable constants with the default coaxial formulation.
"""
function compute(
        problem::CableConstantsProblem;
        options::NamedTuple = (;)
)
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
        options::NamedTuple = (;)
)
    return compute(LineCableModelsCoaxial(), problem, formulation; options)
end

"""
$(TYPEDSIGNATURES)

Compute earth-free cable constants through the coaxial backend tag.
"""
function compute(
        ::LineCableModelsCoaxial,
        problem::CableConstantsProblem,
        formulation::CableConstantsFormulation = CableConstantsFormulation();
        options::NamedTuple = (;)
)
    computation_options(Val(CableConstantsProblem), options)
    validate(problem)
    workspace = CableConstantsWorkspace(problem, formulation)
    return _solve!(workspace, problem, formulation)
end

"""
$(TYPEDSIGNATURES)

Calculate earth-free constants directly from one completed cable design.
"""
function CableConstants(
        design::CableDesign;
        temperature::Real = 20,
        frequency::Real = 1e-3,
        formulation::CableConstantsFormulation = CableConstantsFormulation(),
        options::NamedTuple = (;)
)
    problem = CableConstantsProblem(design; temperature, frequency)
    return compute(problem, formulation; options)
end

computation_owner(::CableConstantsFormulation) = LineCableModelsCoaxial

function computation_details(
        ::Val{LineCableModelsCoaxial},
        ::CableConstants
)::ComputationDetails
    return (;)
end
