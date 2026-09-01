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
end

Base.eltype(::CableConstantsProblem{T}) where {T} = T
Base.eltype(::Type{<:CableConstantsProblem{T}}) where {T} = T

function _check_cable_constants_problem(problem::CableConstantsProblem)
    validate(problem.design)
    isfinite(problem.temperature) || throw(DomainError(
        problem.temperature,
        "cable-constant temperature must be finite",
    ))
    problem.frequency in (oftype(problem.frequency, 50), oftype(problem.frequency, 60)) ||
        throw(DomainError(
            problem.frequency,
            "cable-constant base frequency must be 50 Hz or 60 Hz",
        ))
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
    problem = CableConstantsProblem{T, typeof(design)}(
        design,
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
        internal_impedance = formula(:default),
        insulation_impedance = formula(:default),
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
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

"""
$(TYPEDEF)

Own the local cable arrays, corrected resistivities, and reusable matrices for
one cable-constant calculation. The constructor consumes the structural
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
    cable = LocalCableData(blueprint)
    @inbounds for assembly in cable.assemblies
        isempty(cable.dielectric_ranges[first(assembly)]) && throw(ArgumentError(
            "assembly core :$(cable.terminals[first(assembly)]) has no radial dielectric path",
        ))
    end
    count = length(cable.terminals)
    rho = copy(cable.rho0_cond)
    @inbounds for index in eachindex(rho)
        if formulation.options.temperature_correction
            rho[index] *= one(T) +
                          cable.alpha_cond[index] *
                          (problem.temperature - cable.T0_cond[index])
        end
    end

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
        R = Vector{T}(undef, length(cable.assemblies)),
        L = Vector{T}(undef, length(cable.assemblies)),
        C = Vector{T}(undef, length(cable.assemblies)),
        G = Vector{T}(undef, length(cable.assemblies))
    )
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
    cable_impedance!(
        buffers.Z,
        workspace.cable,
        workspace.rho,
        formulation.methods,
        s
    )
    cable_admittance!(
        buffers.Y,
        workspace.cable,
        formulation.methods,
        problem.frequency,
        problem.temperature,
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
        engine::LineCableModelsCoaxial,
        problem::CableConstantsProblem,
        formulation::CableConstantsFormulation = CableConstantsFormulation();
        options::NamedTuple = (;)
)
    computation_options(Val(CableConstantsProblem), options)
    validate(problem)
    blueprint = flatten(engine, problem.design, eltype(problem))
    workspace = CableConstantsWorkspace(problem, formulation, blueprint)
    return _solve!(workspace, problem, formulation)
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
