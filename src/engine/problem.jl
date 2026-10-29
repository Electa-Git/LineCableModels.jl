"""
$(TYPEDEF)

Describe a line-parameters computation by its cable system, earth, temperature, and
frequencies. The computation converts numerical inputs to `T`, which is at least `Float64`.

$(TYPEDFIELDS)
"""
struct LineParametersProblem{T <: Real, S <: CableSystem, E <: EarthModel{T}}
    "Cables of the line in one cross-section, with the conductor index of each terminal."
    system::S
    "Air and the earth below the cables."
    earth::E
    "Temperature of every conductor and dielectric \\[°C\\]."
    temperature::T
    "Frequencies of the scan in strict order from low to high \\[Hz\\]."
    frequencies::Vector{T}
    function LineParametersProblem{T}(system::S, earth::E, temperature::T,
            frequencies::Vector{T}) where {T <: Real, S <: CableSystem, E <: EarthModel{T}}
        return validate(new{T, S, E}(system, earth, temperature, frequencies))
    end
end

"""
$(TYPEDSIGNATURES)

Build a problem for the cables of `system` above or in `earth`. Promote numerical inputs
with `Float64` to select the computation type.

# Arguments

- `system`: cables and their connections in a cross-section, with the line length.
- `earth`: air and earth material properties.

# Keywords

- `temperature`: temperature of every conductor and dielectric \\[°C\\]. Default: `20`.
- `frequencies`: positive, finite frequencies in strict order from low to high \\[Hz\\].

# Returns

- A [`LineParametersProblem`](@ref).

# Errors

- Throws `DomainError` for an empty frequency vector or a nonfinite temperature.
- Throws `DomainError` for a nonpositive or nonfinite frequency or frequencies outside
  strict order from low to high.
- Throws `ArgumentError` when no terminal connects to a conductor index `k ≥ 1`.
- Throws `ArgumentError` when a cable crosses the interface between air and earth. A cable
  can touch the interface.
"""
function LineParametersProblem(system::CableSystem, earth::EarthModel;
        temperature::Real = 20, frequencies::AbstractVector{<:Real})
    T = promote_type(Float64, eltype(system), eltype(earth), typeof(float(temperature)),
        float(eltype(frequencies)))
    return LineParametersProblem{T}(system, convert(EarthModel{T}, earth),
        convert(T, temperature), convert(Vector{T}, frequencies))
end

# The frequencies increase, the temperature is finite, and a terminal connects to a live
# conductor. A cable at height `y` with outer radius `R` crosses the interface between the
# air and the earth when |y| < R.
function validate(problem::LineParametersProblem)
    frequencies = problem.frequencies
    isempty(frequencies) && throw(DomainError(frequencies,
        "LineParametersProblem.frequencies must not be empty"))
    all(f -> isfinite(f) && f > zero(f), frequencies) || throw(DomainError(frequencies,
        "LineParametersProblem.frequencies must be positive and finite"))
    issorted(frequencies; lt = <=) || throw(DomainError(frequencies,
        "LineParametersProblem.frequencies must increase"))
    isfinite(problem.temperature) || throw(DomainError(problem.temperature,
        "LineParametersProblem.temperature must be finite"))
    system = problem.system
    any(row -> row.connection >= 1, system.terminals) || throw(ArgumentError(
        "the system $(system.id) needs a terminal connected to a conductor k ≥ 1"))
    for (k, design) in enumerate(system.designs)
        y, R = system.placements[k].y, support(design.outer)
        abs(y) < R && throw(ArgumentError(
            "the cable $k of $(system.id) crosses the interface between the air and the " *
            "earth: its center lies $(abs(y)) m from it, closer than its outer radius $R m"))
    end
    return problem
end

function Base.map(f, problem::LineParametersProblem)
    return LineParametersProblem(map(f, problem.system), map(f, problem.earth);
        temperature = f(problem.temperature), frequencies = map(f, problem.frequencies))
end
