"""
    ParametricSweep{C,R} <: AbstractVector{R}

Store an ordered deterministic parameter study as paired cases and results.
Indexing and iteration return results; [`cases`](@ref) preserves the input that
produced each result. The container records no Cartesian axes or GridSpace
state.

All cases must have one concrete type and all results must have one concrete
type. This keeps scalar indexing type-stable while allowing future GridSpace
code to supply the ordered cases.
"""
struct ParametricSweep{C, R} <: AbstractVector{R}
    "Ordered inputs evaluated by the parameter study."
    cases::Vector{C}
    "Ordered results paired one-to-one with `cases`."
    results::Vector{R}

    function ParametricSweep{C, R}(
            cases::Vector{C},
            results::Vector{R}
    ) where {C, R}
        length(cases) == length(results) || throw(
            DimensionMismatch("parametric cases and results must have equal lengths"),
        )
        return new{C, R}(cases, results)
    end
end

function _sweep_vector(name::AbstractString, values)
    collected = collect(values)
    if isempty(collected)
        isconcretetype(eltype(collected)) || throw(
            ArgumentError("empty parametric $name must declare a concrete element type"),
        )
        return collected
    end
    concrete_type = typeof(first(collected))
    all(value -> typeof(value) === concrete_type, collected) || throw(
        ArgumentError("parametric $name must have one concrete type"),
    )
    return Vector{concrete_type}(collected)
end

@doc """
    ParametricSweep(cases, results)

Construct an ordered deterministic parameter study.

# Arguments

- `cases`: Inputs supplied to the deterministic calculations.
- `results`: Results in the same order as `cases`.

# Returns

- A `ParametricSweep` whose scalar indexing and iteration return results.

# Errors

Throws `ArgumentError` when either collection is heterogeneous, or when an
empty collection does not declare a concrete element type.
Throws `DimensionMismatch` when the collections have different lengths.

# Examples

```julia
sweep = ParametricSweep([:case_a, :case_b], [result_a, result_b])
first(sweep) == result_a
cases(sweep)[2] == :case_b
```
"""
function ParametricSweep(cases, results)
    case_values = _sweep_vector("cases", cases)
    result_values = _sweep_vector("results", results)
    return ParametricSweep{eltype(case_values), eltype(result_values)}(
        case_values,
        result_values
    )
end

Base.IndexStyle(::Type{<:ParametricSweep}) = IndexLinear()
Base.size(sweep::ParametricSweep) = (length(sweep.results),)
Base.getindex(sweep::ParametricSweep, index::Int) = sweep.results[index]
function Base.getindex(sweep::ParametricSweep, indices::AbstractVector{<:Integer})
    return ParametricSweep(sweep.cases[indices], sweep.results[indices])
end
function Base.getindex(sweep::ParametricSweep, indices::AbstractRange{<:Integer})
    return ParametricSweep(sweep.cases[indices], sweep.results[indices])
end
function Base.getindex(sweep::ParametricSweep, ::Colon)
    ParametricSweep(copy(sweep.cases), copy(sweep.results))
end
Base.copy(sweep::ParametricSweep) = sweep[:]

"""
    cases(sweep)

Return the ordered inputs paired with the results in `sweep`.
"""
cases(sweep::ParametricSweep) = sweep.cases

"""
    results(sweep)

Return the ordered result vector stored by `sweep`.
"""
results(sweep::ParametricSweep) = sweep.results

"""
    ncases(sweep)

Return the number of paired cases and results in `sweep`.
"""
ncases(sweep::ParametricSweep) = length(sweep)

function _common_result_value(accessor, sweep::ParametricSweep, name::Symbol)
    isempty(sweep) && throw(ArgumentError("$name is unavailable for an empty sweep"))
    value = accessor(first(sweep))
    all(index -> isequal(accessor(sweep[index]), value), 2:length(sweep)) || throw(
        ArgumentError("$name varies between cases; use $name(sweep, case_index)"),
    )
    return value
end

"""
    basis(sweep[, case_index])
    domain(sweep[, case_index])
    frequencies(sweep[, case_index])
    nconductors(sweep[, case_index])
    nfrequencies(sweep[, case_index])

Return a shared result property or the property for one case. The unindexed
form throws when the sweep is empty or the property differs between cases.
"""
basis(sweep::ParametricSweep) = _common_result_value(basis, sweep, :basis)
domain(sweep::ParametricSweep) = _common_result_value(domain, sweep, :domain)
frequencies(sweep::ParametricSweep) = _common_result_value(frequencies, sweep, :frequencies)
nconductors(sweep::ParametricSweep) = _common_result_value(nconductors, sweep, :nconductors)
function nfrequencies(sweep::ParametricSweep)
    _common_result_value(nfrequencies, sweep, :nfrequencies)
end

basis(sweep::ParametricSweep, case_index::Integer) = basis(sweep[Int(case_index)])
domain(sweep::ParametricSweep, case_index::Integer) = domain(sweep[Int(case_index)])
function frequencies(sweep::ParametricSweep, case_index::Integer)
    frequencies(sweep[Int(case_index)])
end
function nconductors(sweep::ParametricSweep, case_index::Integer)
    nconductors(sweep[Int(case_index)])
end
function nfrequencies(sweep::ParametricSweep, case_index::Integer)
    nfrequencies(sweep[Int(case_index)])
end

"""
    Z(sweep[, case_index[, i, j[, k]]])
    Y(sweep[, case_index[, i, j[, k]]])
    R(sweep[, case_index[, i, j[, k]]])
    X(sweep[, case_index[, i, j[, k]]])
    L(sweep[, case_index[, i, j[, k]]])
    G(sweep[, case_index[, i, j[, k]]])
    B(sweep[, case_index[, i, j[, k]]])
    C(sweep[, case_index[, i, j[, k]]])

Return a component for every result or for one selected case. After
`case_index`, selection follows the accessor grammar of the stored result. For
`LineParameters`, `(i, j)` returns the complete frequency response and `k` may
be a scalar, range, or `:`. Canonical units follow [`basis`](@ref).
"""
Z(sweep::ParametricSweep) = map(Z, sweep.results)
Y(sweep::ParametricSweep) = map(Y, sweep.results)
R(sweep::ParametricSweep) = map(R, sweep.results)
X(sweep::ParametricSweep) = map(X, sweep.results)
L(sweep::ParametricSweep) = map(L, sweep.results)
G(sweep::ParametricSweep) = map(G, sweep.results)
B(sweep::ParametricSweep) = map(B, sweep.results)
C(sweep::ParametricSweep) = map(C, sweep.results)

Z(sweep::ParametricSweep, case_index::Integer, args...) = Z(sweep[Int(case_index)], args...)
Y(sweep::ParametricSweep, case_index::Integer, args...) = Y(sweep[Int(case_index)], args...)
R(sweep::ParametricSweep, case_index::Integer, args...) = R(sweep[Int(case_index)], args...)
X(sweep::ParametricSweep, case_index::Integer, args...) = X(sweep[Int(case_index)], args...)
L(sweep::ParametricSweep, case_index::Integer, args...) = L(sweep[Int(case_index)], args...)
G(sweep::ParametricSweep, case_index::Integer, args...) = G(sweep[Int(case_index)], args...)
B(sweep::ParametricSweep, case_index::Integer, args...) = B(sweep[Int(case_index)], args...)
C(sweep::ParametricSweep, case_index::Integer, args...) = C(sweep[Int(case_index)], args...)

series_impedance(sweep::ParametricSweep) = map(series_impedance, sweep.results)
shunt_admittance(sweep::ParametricSweep) = map(shunt_admittance, sweep.results)
function series_impedance(sweep::ParametricSweep, case_index::Integer)
    series_impedance(sweep[Int(case_index)])
end
function shunt_admittance(sweep::ParametricSweep, case_index::Integer)
    shunt_admittance(sweep[Int(case_index)])
end
resistance(sweep::ParametricSweep, args...) = R(sweep, args...)
reactance(sweep::ParametricSweep, args...) = X(sweep, args...)
inductance(sweep::ParametricSweep, args...) = L(sweep, args...)
conductance(sweep::ParametricSweep, args...) = G(sweep, args...)
susceptance(sweep::ParametricSweep, args...) = B(sweep, args...)
capacitance(sweep::ParametricSweep, args...) = C(sweep, args...)
