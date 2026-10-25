"""
$(TYPEDEF)

Describe one computation of line parameters: the cables of the line, the earth below them,
the operating temperature and the frequencies of the scan. The problem keeps the system and
the earth as given. `T` is the floating type of the temperature and the frequencies, and the
computation converts the values that it reads to `T`.

$(TYPEDFIELDS)
"""
struct LineParametersProblem{T <: Real, S <: CableSystem, E <: EarthModel}
    "Cables of the line in one cross-section, with the conductor index of each terminal."
    system::S
    "Air and the earth below the cables."
    earth::E
    "Operating temperature of every conductor and dielectric \\[°C\\]."
    temperature::T
    "Frequencies of the scan, increasing \\[Hz\\]."
    frequencies::Vector{T}
    function LineParametersProblem{T}(system::S, earth::E, temperature::T,
            frequencies::Vector{T}) where {T <: Real, S <: CableSystem, E <: EarthModel}
        return validate(new{T, S, E}(system, earth, temperature, frequencies))
    end
end

"""
$(TYPEDSIGNATURES)

Build the problem of the cables of `system` above or in `earth`. The temperature and the
frequencies take the floating type that the system, the earth, the temperature and the
frequencies promote to.

# Arguments

- `system`: cables of the line in one cross-section.
- `earth`: air and the earth below the cables.

# Keywords

- `temperature`: operating temperature of every conductor and dielectric \\[°C\\]. Default:
  `20`.
- `frequencies`: frequencies of the scan, positive, finite and increasing \\[Hz\\].

# Returns

- A [`LineParametersProblem`](@ref).

# Errors

- Throws `DomainError` when the frequencies are empty, not positive and finite, or not
  increasing, or when the temperature is not finite.
- Throws `ArgumentError` when no terminal connects to a conductor `k ≥ 1`, or when the
  outline of a cable crosses the interface between the air and the earth. A cable can
  touch the interface.
"""
function LineParametersProblem(system::CableSystem, earth::EarthModel;
        temperature::Real = 20, frequencies::AbstractVector{<:Real})
    T = promote_type(eltype(system), eltype(earth), typeof(float(temperature)),
        float(eltype(frequencies)))
    return LineParametersProblem{T}(
        system, earth, convert(T, temperature), convert(Vector{T}, frequencies))
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
