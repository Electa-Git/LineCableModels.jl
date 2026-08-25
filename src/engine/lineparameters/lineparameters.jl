const LINE_PARAMETER_BASES = (:pul, :total)

@inline function _check_basis(value::Symbol)
    value in LINE_PARAMETER_BASES || throw(
        ArgumentError("basis must be :pul or :total; got :$value"),
    )
    return value
end

"""
    SeriesImpedance{T, Basis}

Store a square series-impedance matrix over frequency. Values use \\[Ω/m\\]
when `Basis` is `:pul` and \\[Ω\\] when it is `:total`.
"""
struct SeriesImpedance{T, Basis} <: AbstractArray{T, 3}
    "Complex series-impedance tensor with dimensions conductor × conductor × frequency."
    values::Array{T, 3}

    function SeriesImpedance{T, Basis}(values::Array{T, 3}) where {T, Basis}
        Basis isa Symbol || throw(
            ArgumentError("basis must be :pul or :total; got $(repr(Basis))"),
        )
        _check_basis(Basis)
        return new{T, Basis}(values)
    end
end

"""
    ShuntAdmittance{T, Basis}

Store a square shunt-admittance matrix over frequency. Values use \\[S/m\\]
when `Basis` is `:pul` and \\[S\\] when it is `:total`.
"""
struct ShuntAdmittance{T, Basis} <: AbstractArray{T, 3}
    "Complex shunt-admittance tensor with dimensions conductor × conductor × frequency."
    values::Array{T, 3}

    function ShuntAdmittance{T, Basis}(values::Array{T, 3}) where {T, Basis}
        Basis isa Symbol || throw(
            ArgumentError("basis must be :pul or :total; got $(repr(Basis))"),
        )
        _check_basis(Basis)
        return new{T, Basis}(values)
    end
end

function SeriesImpedance(A::AbstractArray{T, 3}; basis::Symbol = :pul) where {T}
    _check_basis(basis)
    return SeriesImpedance{T, basis}(Array(A))
end

function ShuntAdmittance(A::AbstractArray{T, 3}; basis::Symbol = :pul) where {T}
    _check_basis(basis)
    return ShuntAdmittance{T, basis}(Array(A))
end

@inline basis(::Type{<:SeriesImpedance{T, Basis}}) where {T, Basis} = Basis
@inline basis(::SeriesImpedance{T, Basis}) where {T, Basis} = Basis
@inline basis(::Type{<:ShuntAdmittance{T, Basis}}) where {T, Basis} = Basis
@inline basis(::ShuntAdmittance{T, Basis}) where {T, Basis} = Basis

"""
    LineParameters{T, U, D, Basis}

Frequency-dependent series-impedance and shunt-admittance matrices.

`Basis` is either `:pul` or `:total`. Per-length values are stored in
Ω/m and S/m. Total values are stored in Ω and S. Frequencies are stored in Hz.
"""
struct LineParameters{
    T <: Complex,
    U <: Real,
    D <: LineParamsDomain,
    Basis
} <: AbstractProblemResult
    "Frequency-dependent series impedance \\[Ω/m\\] or \\[Ω\\]."
    Z::SeriesImpedance{T, Basis}
    "Frequency-dependent shunt admittance \\[S/m\\] or \\[S\\]."
    Y::ShuntAdmittance{T, Basis}
    "Frequency samples \\[Hz\\]."
    f::Vector{U}

    function LineParameters(
            ::Type{D},
            Z::SeriesImpedance{T, Basis},
            Y::ShuntAdmittance{T, Basis},
            f::AbstractVector{U}
    ) where {
            D <: LineParamsDomain,
            T <: Complex,
            U <: Real,
            Basis
    }
        _check_basis(Basis)
        size(Z, 1) == size(Z, 2) || throw(DimensionMismatch("Z must be square"))
        size(Y, 1) == size(Y, 2) || throw(DimensionMismatch("Y must be square"))
        size(Z) == size(Y) || throw(
            DimensionMismatch("Z and Y must have equal n×n×nfreq dimensions"),
        )
        size(Z, 3) == length(f) || throw(
            DimensionMismatch("frequency count must match the Z/Y third dimension"),
        )
        all(isfinite, f) || throw(ArgumentError("frequencies must be finite"))
        return new{T, U, D, Basis}(Z, Y, Vector{U}(f))
    end
end

function LineParameters(
        Z::SeriesImpedance{T, Basis},
        Y::ShuntAdmittance{T, Basis},
        f::AbstractVector{U}
) where {T <: Complex, U <: Real, Basis}
    return LineParameters(PhaseDomain, Z, Y, f)
end

function LineParameters(
        Z::SeriesImpedance{TZ, ZBasis},
        Y::ShuntAdmittance{TY, YBasis},
        f::AbstractVector{U}
) where {
        TZ <: Complex,
        TY <: Complex,
        U <: Real,
        ZBasis,
        YBasis
}
    ZBasis === YBasis || throw(
        ArgumentError("Z and Y must have the same basis; got :$ZBasis and :$YBasis"),
    )
    element_type = promote_type(TZ, TY)
    return LineParameters(
        PhaseDomain,
        SeriesImpedance(convert(Array{element_type, 3}, Z.values); basis = ZBasis),
        ShuntAdmittance(convert(Array{element_type, 3}, Y.values); basis = YBasis),
        f
    )
end

function LineParameters(
        ::Type{D},
        Z::AbstractArray{TZ, 3},
        Y::AbstractArray{TY, 3},
        f::AbstractVector{U};
        basis::Symbol = :pul
) where {
        D <: LineParamsDomain,
        TZ <: Complex,
        TY <: Complex,
        U <: Real
}
    _check_basis(basis)
    element_type = promote_type(TZ, TY)
    return LineParameters(
        D,
        SeriesImpedance(convert(Array{element_type, 3}, Z); basis),
        ShuntAdmittance(convert(Array{element_type, 3}, Y); basis),
        f
    )
end

function LineParameters(
        Z::AbstractArray{TZ, 3},
        Y::AbstractArray{TY, 3},
        f::AbstractVector{U};
        basis::Symbol = :pul
) where {
        TZ <: Complex,
        TY <: Complex,
        U <: Real
}
    return LineParameters(PhaseDomain, Z, Y, f; basis)
end

@inline domain(::Type{<:LineParameters{T, U, D}}) where {T, U, D <: LineParamsDomain} = D
@inline domain(lp::LineParameters) = domain(typeof(lp))
@inline basis(::Type{<:LineParameters{T, U, D, Basis}}) where {T, U, D, Basis} = Basis
@inline basis(::LineParameters{T, U, D, Basis}) where {T, U, D, Basis} = Basis

observe(lp::LineParameters, ::typeof(frequencies)) = lp.f
observe(lp::LineParameters, ::typeof(frequencies), indices...) = getindex(lp.f, indices...)
frequencies(lp::LineParameters, indices...) = observe(lp, frequencies, indices...)
nconductors(lp::LineParameters) = size(lp.Z, 1)
nfrequencies(lp::LineParameters) = length(lp.f)

observables(::Type{<:LineParameters}) = (
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
    (Y, angle)
)

series_impedance(value::Union{LineParameters, SeriesImpedance}) = Z(value)
shunt_admittance(value::Union{LineParameters, ShuntAdmittance}) = Y(value)
"""
    Z(parameters[, i, j[, k]])
    Y(parameters[, i, j[, k]])

Return series impedance or shunt admittance. With `(i, j)`, return the complete
frequency response at that matrix position. `k` may be one index, a range, or
`:`. Stored units follow [`basis`](@ref): \\[Ω/m\\] and \\[S/m\\] for
`:pul`, or \\[Ω\\] and \\[S\\] for `:total`.
"""
@inline _observe_array(values::AbstractArray) = values
@inline _observe_array(values::AbstractArray, i, j) = view(values, i, j, :)
@inline _observe_array(values::AbstractArray, indices...) = getindex(values, indices...)

observe(impedance::SeriesImpedance, ::typeof(Z), indices...) =
    _observe_array(impedance.values, indices...)
observe(admittance::ShuntAdmittance, ::typeof(Y), indices...) =
    _observe_array(admittance.values, indices...)
observe(lp::LineParameters, ::typeof(Z), indices...) =
    _observe_array(lp.Z.values, indices...)
observe(lp::LineParameters, ::typeof(Y), indices...) =
    _observe_array(lp.Y.values, indices...)

observe(value::Union{LineParameters, SeriesImpedance}, ::typeof(R), indices...) =
    real.(observe(value, Z, indices...))
observe(value::Union{LineParameters, SeriesImpedance}, ::typeof(X), indices...) =
    imag.(observe(value, Z, indices...))
observe(value::Union{LineParameters, ShuntAdmittance}, ::typeof(G), indices...) =
    real.(observe(value, Y, indices...))
observe(value::Union{LineParameters, ShuntAdmittance}, ::typeof(B), indices...) =
    imag.(observe(value, Y, indices...))

observe(value::LineParameters, ::typeof(Z), ::typeof(abs), indices...) =
    abs.(observe(value, Z, indices...))
observe(value::SeriesImpedance, ::typeof(Z), ::typeof(abs), indices...) =
    abs.(observe(value, Z, indices...))
observe(value::LineParameters, ::typeof(Z), ::typeof(angle), indices...) =
    angle.(observe(value, Z, indices...))
observe(value::SeriesImpedance, ::typeof(Z), ::typeof(angle), indices...) =
    angle.(observe(value, Z, indices...))
observe(value::LineParameters, ::typeof(Y), ::typeof(abs), indices...) =
    abs.(observe(value, Y, indices...))
observe(value::ShuntAdmittance, ::typeof(Y), ::typeof(abs), indices...) =
    abs.(observe(value, Y, indices...))
observe(value::LineParameters, ::typeof(Y), ::typeof(angle), indices...) =
    angle.(observe(value, Y, indices...))
observe(value::ShuntAdmittance, ::typeof(Y), ::typeof(angle), indices...) =
    angle.(observe(value, Y, indices...))

observables(::Type{<:SeriesImpedance}) = (Z, R, X, L, (Z, abs), (Z, angle))
observables(::Type{<:ShuntAdmittance}) = (Y, G, B, C, (Y, abs), (Y, angle))

Z(value::Union{LineParameters, SeriesImpedance}, indices...) = observe(value, Z, indices...)
Y(value::Union{LineParameters, ShuntAdmittance}, indices...) = observe(value, Y, indices...)
R(value::Union{LineParameters, SeriesImpedance}, indices...) = observe(value, R, indices...)
X(value::Union{LineParameters, SeriesImpedance}, indices...) = observe(value, X, indices...)
G(value::Union{LineParameters, ShuntAdmittance}, indices...) = observe(value, G, indices...)
B(value::Union{LineParameters, ShuntAdmittance}, indices...) = observe(value, B, indices...)

resistance(value::Union{LineParameters, SeriesImpedance}, args...) = R(value, args...)
reactance(value::Union{LineParameters, SeriesImpedance}, args...) = X(value, args...)
conductance(value::Union{LineParameters, ShuntAdmittance}, args...) = G(value, args...)
susceptance(value::Union{LineParameters, ShuntAdmittance}, args...) = B(value, args...)

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
        admittance::ShuntAdmittance,
        ::typeof(C),
        frequencies::AbstractVector
)
    size(admittance, 3) == length(frequencies) || throw(
        DimensionMismatch("frequency count must match the admittance third dimension"),
    )
    return imag.(admittance.values) ./ _angular_frequencies(frequencies)
end

"""
$(TYPEDSIGNATURES)

Return series inductance using the same frequency-selection grammar as
[`Z`](@ref), evaluated as:

```math
L(f) = \\frac{\\operatorname{Im} Z(f)}{2\\pi f}.
```

Units are \\[H/m\\] for `:pul` and \\[H\\] for `:total`.

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

Units are \\[F/m\\] for `:pul` and \\[F\\] for `:total`.

# Errors

Throws `DomainError` when a selected frequency is zero.
"""
function observe(lp::LineParameters, ::typeof(C))
    any(iszero, lp.f) && throw(DomainError(lp.f, "C is undefined at zero frequency"))
    return imag.(lp.Y.values) ./ reshape(2π .* lp.f, 1, 1, :)
end

observe(lp::LineParameters, ::typeof(L), i, j) = observe(lp, L, i, j, :)
observe(lp::LineParameters, ::typeof(C), i, j) = observe(lp, C, i, j, :)
observe(lp::LineParameters, ::typeof(L), i, j, k) =
    observe(lp, X, i, j, k) ./ _angular_frequencies(lp, k)
observe(lp::LineParameters, ::typeof(C), i, j, k) =
    observe(lp, B, i, j, k) ./ _angular_frequencies(lp, k)

L(lp::LineParameters, args...) = observe(lp, L, args...)
C(lp::LineParameters, args...) = observe(lp, C, args...)
inductance(lp::LineParameters, args...) = L(lp, args...)
capacitance(lp::LineParameters, args...) = C(lp, args...)
