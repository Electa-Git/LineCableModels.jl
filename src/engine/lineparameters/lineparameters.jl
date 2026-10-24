"""
    SeriesImpedance{T}

Store a square series-impedance matrix over frequency, in \\[Ω/m\\].
"""
struct SeriesImpedance{T} <: AbstractArray{T, 3}
    "Complex series-impedance tensor with dimensions conductor × conductor × frequency."
    values::Array{T, 3}

    SeriesImpedance{T}(values::Array{T, 3}) where {T} = new{T}(values)
end

"""
    ShuntAdmittance{T}

Store a square shunt-admittance matrix over frequency, in \\[S/m\\].
"""
struct ShuntAdmittance{T} <: AbstractArray{T, 3}
    "Complex shunt-admittance tensor with dimensions conductor × conductor × frequency."
    values::Array{T, 3}

    ShuntAdmittance{T}(values::Array{T, 3}) where {T} = new{T}(values)
end

SeriesImpedance(A::AbstractArray{T, 3}) where {T} = SeriesImpedance{T}(Array(A))
ShuntAdmittance(A::AbstractArray{T, 3}) where {T} = ShuntAdmittance{T}(Array(A))

"""
    LineParameters{T, U, D}

Frequency-dependent series-impedance and shunt-admittance matrices of a line, per meter,
with Z in Ω/m and Y in S/m. Frequencies are stored in Hz. `D` is the concrete physical
coordinate-domain value. The record of the computation, [`ComputationDetails`](@ref), has
one type for every content. The whole line's matrices are `Z(result) .* line_length(result)`
and `Y(result) .* line_length(result)`.
"""
struct LineParameters{
    T <: Complex,
    U <: Real,
    D <: LineParamsDomain
} <: AbstractCoreResult
    "Frequency-dependent series impedance \\[Ω/m\\]."
    Z::SeriesImpedance{T}
    "Frequency-dependent shunt admittance \\[S/m\\]."
    Y::ShuntAdmittance{T}
    "Frequency samples \\[Hz\\]."
    f::Vector{U}
    "Physical coordinate domain, including any required transformation state."
    domain::D
    "The line's terminals after reduction, one per matrix row of the phase domain."
    terminals::Vector{String}
    "Length of the line \\[m\\]."
    line_length::U
    "Supplemental output retained by the computation."
    details::ComputationDetails

    function LineParameters(
            domain::D,
            Z::SeriesImpedance{T},
            Y::ShuntAdmittance{T},
            f::AbstractVector{U},
            terminals::AbstractVector{<:AbstractString},
            line_length::Real,
            details::ComputationDetails
    ) where {
            D <: LineParamsDomain,
            T <: Complex,
            U <: Real
    }
        retained=haskey(details.data,:gridpoint) ? details :
            ComputationDetails(merge(details.data,(gridpoint=Commons.gridpoint_id(),)))
        return validate(new{T, U, D}(
            Z,
            Y,
            Vector{U}(f),
            domain,
            terminals,
            line_length,
            retained
        ))
    end
end

# The phase domain does not restrict the coefficients.
validate(domain::LineParamsDomain, ::LineParameters) = domain

function validate(domain::ModalDomain, parameters::LineParameters)
    size(domain.operators)==size(parameters.Z.values) || throw(DimensionMismatch(
        "modal operators must align with line-parameter coefficients"))
    return domain
end

function validate(parameters::LineParameters)
    size(parameters.Z, 1) == size(parameters.Z, 2) || throw(DimensionMismatch(
        "LineParameters.Z must be square; received size $(size(parameters.Z))"
    ))
    size(parameters.Y, 1) == size(parameters.Y, 2) || throw(DimensionMismatch(
        "LineParameters.Y must be square; received size $(size(parameters.Y))"
    ))
    size(parameters.Z) == size(parameters.Y) || throw(DimensionMismatch(
        "LineParameters.Z and LineParameters.Y must have equal n×n×nfreq " *
        "dimensions; received $(size(parameters.Z)) and $(size(parameters.Y))"
    ))
    size(parameters.Z, 3) == length(parameters.f) || throw(DimensionMismatch(
        "LineParameters.f must contain one value per matrix frequency plane; " *
        "received $(length(parameters.f)) values for $(size(parameters.Z, 3)) planes"
    ))
    validate(parameters.domain, parameters)
    all(isfinite, parameters.f) || throw(ArgumentError(
        "LineParameters.f must contain only finite frequencies; received " *
        repr(parameters.f)
    ))
    length(parameters.terminals) == size(parameters.Z, 1) || throw(DimensionMismatch(
        "LineParameters.terminals must name one terminal per matrix row; received " *
        "$(length(parameters.terminals)) names for $(size(parameters.Z, 1)) rows"
    ))
    isfinite(parameters.line_length) && parameters.line_length > zero(parameters.line_length) ||
        throw(DomainError(parameters.line_length,
            "LineParameters.line_length must be positive and finite"))
    return parameters
end

# Without `terminals`, a result built from matrices names its terminals by position. Its
# line length is always given.
function LineParameters(
        domain::LineParamsDomain,
        Z::SeriesImpedance{T},
        Y::ShuntAdmittance{T},
        f::AbstractVector{<:Real};
        terminals::AbstractVector{<:AbstractString} = string.(1:size(Z, 1)),
        line_length::Real,
        details::ComputationDetails = ComputationDetails()
) where {T <: Complex}
    return LineParameters(domain, Z, Y, f, terminals, line_length, details)
end

function LineParameters(
        ::Type{PhaseDomain},
        Z::SeriesImpedance{T},
        Y::ShuntAdmittance{T},
        f::AbstractVector{<:Real};
        kwargs...
) where {T <: Complex}
    return LineParameters(PhaseDomain(), Z, Y, f; kwargs...)
end

function LineParameters(
        Z::SeriesImpedance{T},
        Y::ShuntAdmittance{T},
        f::AbstractVector{<:Real};
        kwargs...
) where {T <: Complex}
    return LineParameters(PhaseDomain, Z, Y, f; kwargs...)
end

function LineParameters(
        Z::SeriesImpedance{TZ},
        Y::ShuntAdmittance{TY},
        f::AbstractVector{U};
        kwargs...
) where {
        TZ <: Complex,
        TY <: Complex,
        U <: Real
}
    element_type = promote_type(TZ, TY)
    return LineParameters(
        PhaseDomain,
        SeriesImpedance(convert(Array{element_type, 3}, Z.values)),
        ShuntAdmittance(convert(Array{element_type, 3}, Y.values)),
        f;
        kwargs...
    )
end

function LineParameters(
        domain::D,
        Z::AbstractArray{TZ, 3},
        Y::AbstractArray{TY, 3},
        f::AbstractVector{U};
        kwargs...
) where {
        D <: LineParamsDomain,
        TZ <: Complex,
        TY <: Complex,
        U <: Real
}
    element_type = promote_type(TZ, TY)
    return LineParameters(
        domain,
        SeriesImpedance(convert(Array{element_type, 3}, Z)),
        ShuntAdmittance(convert(Array{element_type, 3}, Y)),
        f;
        kwargs...
    )
end

function LineParameters(
        ::Type{PhaseDomain},
        Z::AbstractArray{TZ, 3},
        Y::AbstractArray{TY, 3},
        f::AbstractVector{U};
        kwargs...
) where {TZ <: Complex, TY <: Complex, U <: Real}
    return LineParameters(PhaseDomain(), Z, Y, f; kwargs...)
end

function LineParameters(
        Z::AbstractArray{TZ, 3},
        Y::AbstractArray{TY, 3},
        f::AbstractVector{U};
        kwargs...
) where {
        TZ <: Complex,
        TY <: Complex,
        U <: Real
}
    return LineParameters(PhaseDomain, Z, Y, f; kwargs...)
end

@inline domain(::Type{<:LineParameters{
    T, U, D}}) where {T, U, D <: LineParamsDomain} = domain(D)
@inline domain(lp::LineParameters) = domain(typeof(lp))
details(parameters::LineParameters) = parameters.details
line_length(parameters::LineParameters) = parameters.line_length

observe(lp::LineParameters, ::typeof(frequencies)) = lp.f
observe(lp::LineParameters, ::typeof(frequencies), indices...) = getindex(lp.f, indices...)
frequencies(lp::LineParameters, indices...) = observe(lp, frequencies, indices...)
nconductors(lp::LineParameters) = size(lp.Z, 1)
nfrequencies(lp::LineParameters) = length(lp.f)

function observables(::Type{<:LineParameters})
    (
        frequencies,
        Z,
        Y,
        R,
        X,
        L,
        G,
        B,
        C,
        (Z, abs),
        (Z, angle),
        (Y, abs),
        (Y, angle),
        (Z, diag),
        (Y, diag),
        (R, diag),
        (X, diag),
        (L, diag),
        (G, diag),
        (B, diag),
        (C, diag)
    )
end

series_impedance(value::Union{LineParameters, SeriesImpedance}) = observe(value, Z)
shunt_admittance(value::Union{LineParameters, ShuntAdmittance}) = observe(value, Y)
"""
    Z(parameters[, i, j[, k]])
    Y(parameters[, i, j[, k]])

Return series impedance or shunt admittance. With `(i, j)`, return the complete
frequency response at that matrix position. Use one index, a range or `:` for `k`. Units
are \\[Ω/m\\] and \\[S/m\\].
"""
@inline _observe_array(values::AbstractArray) = values
@inline _observe_array(values::AbstractArray, i, j) = view(values, i, j, :)
@inline _observe_array(values::AbstractArray, indices...) = getindex(values, indices...)

function observe(impedance::SeriesImpedance, ::typeof(Z), indices...)
    _observe_array(impedance.values, indices...)
end
function observe(admittance::ShuntAdmittance, ::typeof(Y), indices...)
    _observe_array(admittance.values, indices...)
end
function observe(lp::LineParameters, ::typeof(Z), indices...)
    _observe_array(lp.Z.values, indices...)
end
function observe(lp::LineParameters, ::typeof(Y), indices...)
    _observe_array(lp.Y.values, indices...)
end

function observe(value::Union{LineParameters, SeriesImpedance}, ::typeof(R), indices...)
    real.(observe(value, Z, indices...))
end
function observe(value::Union{LineParameters, SeriesImpedance}, ::typeof(X), indices...)
    imag.(observe(value, Z, indices...))
end
function observe(value::Union{LineParameters, ShuntAdmittance}, ::typeof(G), indices...)
    real.(observe(value, Y, indices...))
end
function observe(value::Union{LineParameters, ShuntAdmittance}, ::typeof(B), indices...)
    imag.(observe(value, Y, indices...))
end

function observe(value::LineParameters, ::typeof(Z), ::typeof(abs), indices...)
    abs.(observe(value, Z, indices...))
end
function observe(value::SeriesImpedance, ::typeof(Z), ::typeof(abs), indices...)
    abs.(observe(value, Z, indices...))
end
function observe(value::LineParameters, ::typeof(Z), ::typeof(angle), indices...)
    angle.(observe(value, Z, indices...))
end
function observe(value::SeriesImpedance, ::typeof(Z), ::typeof(angle), indices...)
    angle.(observe(value, Z, indices...))
end
function observe(value::LineParameters, ::typeof(Y), ::typeof(abs), indices...)
    abs.(observe(value, Y, indices...))
end
function observe(value::ShuntAdmittance, ::typeof(Y), ::typeof(abs), indices...)
    abs.(observe(value, Y, indices...))
end
function observe(value::LineParameters, ::typeof(Y), ::typeof(angle), indices...)
    angle.(observe(value, Y, indices...))
end
function observe(value::ShuntAdmittance, ::typeof(Y), ::typeof(angle), indices...)
    angle.(observe(value, Y, indices...))
end

function _observe_diagonal(values::AbstractArray{T, 3}, indices...) where {T}
    size(values, 1) == size(values, 2) || throw(DimensionMismatch(
        "diagonal observation requires square matrix slices"
    ))
    diagonal = Matrix{T}(undef, size(values, 1), size(values, 3))
    for sample in axes(values, 3), mode in axes(values, 1)

        diagonal[mode, sample] = values[mode, mode, sample]
    end
    return isempty(indices) ? diagonal : getindex(diagonal, indices...)
end

const _SeriesResult = Union{LineParameters, SeriesImpedance}
const _ShuntResult = Union{LineParameters, ShuntAdmittance}

function observe(value::LineParameters, ::typeof(Z), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, Z), indices...)
end
function observe(value::SeriesImpedance, ::typeof(Z), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, Z), indices...)
end
function observe(value::_SeriesResult, ::typeof(R), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, R), indices...)
end
function observe(value::_SeriesResult, ::typeof(X), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, X), indices...)
end
function observe(value::LineParameters, ::typeof(L), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, L), indices...)
end
"""
$(TYPEDSIGNATURES)

Return diagonal inductances from a standalone impedance tensor and its explicit
frequency vector \\[Hz\\]. The output axes are conductor × frequency. Optional
indices select those axes. Values use \\[H/m\\].
"""
function observe(value::SeriesImpedance, ::typeof(L), ::typeof(diag),
        frequencies::AbstractVector, indices...)
    _observe_diagonal(observe(value, L, frequencies), indices...)
end

function observe(value::LineParameters, ::typeof(Y), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, Y), indices...)
end
function observe(value::ShuntAdmittance, ::typeof(Y), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, Y), indices...)
end
function observe(value::_ShuntResult, ::typeof(G), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, G), indices...)
end
function observe(value::_ShuntResult, ::typeof(B), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, B), indices...)
end
function observe(value::LineParameters, ::typeof(C), ::typeof(diag), indices...)
    _observe_diagonal(observe(value, C), indices...)
end
"""
$(TYPEDSIGNATURES)

Return diagonal capacitances from a standalone admittance tensor and its explicit
frequency vector \\[Hz\\]. The output axes are conductor × frequency. Optional
indices select those axes. Values use \\[F/m\\].
"""
function observe(value::ShuntAdmittance, ::typeof(C), ::typeof(diag),
        frequencies::AbstractVector, indices...)
    _observe_diagonal(observe(value, C, frequencies), indices...)
end

function observables(::Type{<:SeriesImpedance})
    (Z, R, X, L, (Z, abs), (Z, angle), (Z, diag), (R, diag), (X, diag), (L, diag))
end
function observables(::Type{<:ShuntAdmittance})
    (Y, G, B, C, (Y, abs), (Y, angle), (Y, diag), (G, diag), (B, diag), (C, diag))
end

Z(value::Union{LineParameters, SeriesImpedance}, indices...) = observe(value, Z, indices...)
Y(value::Union{LineParameters, ShuntAdmittance}, indices...) = observe(value, Y, indices...)
R(value::Union{LineParameters, SeriesImpedance}, indices...) = observe(value, R, indices...)
X(value::Union{LineParameters, SeriesImpedance}, indices...) = observe(value, X, indices...)
G(value::Union{LineParameters, ShuntAdmittance}, indices...) = observe(value, G, indices...)
B(value::Union{LineParameters, ShuntAdmittance}, indices...) = observe(value, B, indices...)

resistance(value::Union{LineParameters, SeriesImpedance}, args...) = observe(value, R, args...)
reactance(value::Union{LineParameters, SeriesImpedance}, args...) = observe(value, X, args...)
conductance(value::Union{LineParameters, ShuntAdmittance}, args...) = observe(value, G, args...)
susceptance(value::Union{LineParameters, ShuntAdmittance}, args...) = observe(value, B, args...)

@inline function _angular_frequencies(lp::LineParameters, k)
    selected = lp.f[k]
    any(iszero, selected isa Number ? (selected,) : selected) && throw(
        DomainError(selected, "L and C are undefined at zero frequency"),
    )
    return 2π .* selected
end

function _angular_frequencies(frequencies::AbstractVector)
    any(iszero, frequencies) && throw(
        DomainError(frequencies, "L and C are undefined at zero frequency"),
    )
    return reshape(2π .* frequencies, 1, 1, :)
end

function observe(
        impedance::SeriesImpedance,
        ::typeof(L),
        frequencies::AbstractVector
)
    size(impedance, 3) == length(frequencies) || throw(
        DimensionMismatch("frequency count must match the impedance third dimension"),
    )
    return imag.(impedance.values) ./ _angular_frequencies(frequencies)
end

function observe(
        impedance::SeriesImpedance,
        ::typeof(L),
        frequencies::AbstractVector,
        indices...
)
    return getindex(observe(impedance, L, frequencies), indices...)
end

function observe(
        admittance::ShuntAdmittance,
        ::typeof(C),
        frequencies::AbstractVector
)
    size(admittance, 3) == length(frequencies) || throw(
        DimensionMismatch("frequency count must match the admittance third dimension"),
    )
    return imag.(admittance.values) ./ _angular_frequencies(frequencies)
end

function observe(
        admittance::ShuntAdmittance,
        ::typeof(C),
        frequencies::AbstractVector,
        indices...
)
    return getindex(observe(admittance, C, frequencies), indices...)
end

"""
$(TYPEDSIGNATURES)

Return series inductance using the same frequency-selection grammar as
[`Z`](@ref), evaluated as:

```math
L(f) = \\frac{\\operatorname{Im} Z(f)}{2\\pi f}.
```

Units are \\[H/m\\].

# Errors

Throws `DomainError` when a selected frequency is zero.
"""
function observe(lp::LineParameters, ::typeof(L))
    any(iszero, lp.f) && throw(DomainError(lp.f, "L is undefined at zero frequency"))
    return imag.(lp.Z.values) ./ reshape(2π .* lp.f, 1, 1, :)
end

"""
$(TYPEDSIGNATURES)

Return shunt capacitance using the same frequency-selection grammar as
[`Y`](@ref), evaluated as:

```math
C(f) = \\frac{\\operatorname{Im} Y(f)}{2\\pi f}.
```

Units are \\[F/m\\].

# Errors

Throws `DomainError` when a selected frequency is zero.
"""
function observe(lp::LineParameters, ::typeof(C))
    any(iszero, lp.f) && throw(DomainError(lp.f, "C is undefined at zero frequency"))
    return imag.(lp.Y.values) ./ reshape(2π .* lp.f, 1, 1, :)
end

observe(lp::LineParameters, ::typeof(L), i, j) = observe(lp, L, i, j, :)
observe(lp::LineParameters, ::typeof(C), i, j) = observe(lp, C, i, j, :)
function _divide_by_angular_frequency(component, angular_frequency)
    angular_frequency isa Number && return component ./ angular_frequency
    ndims(component) == 1 && return component ./ angular_frequency
    dimensions = (ntuple(_ -> 1, ndims(component) - 1)..., length(angular_frequency))
    return component ./ reshape(angular_frequency, dimensions)
end

function observe(lp::LineParameters, ::typeof(L), i, j, k)
    return _divide_by_angular_frequency(
        observe(lp, X, i, j, k),
        _angular_frequencies(lp, k)
    )
end
function observe(lp::LineParameters, ::typeof(C), i, j, k)
    return _divide_by_angular_frequency(
        observe(lp, B, i, j, k),
        _angular_frequencies(lp, k)
    )
end

L(lp::LineParameters, args...) = observe(lp, L, args...)
C(lp::LineParameters, args...) = observe(lp, C, args...)
inductance(lp::LineParameters, args...) = observe(lp, L, args...)
capacitance(lp::LineParameters, args...) = observe(lp, C, args...)
