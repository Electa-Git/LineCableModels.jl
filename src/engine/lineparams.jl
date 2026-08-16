const LINE_PARAMETER_BASES = (:per_length, :total)

@inline function _check_basis(value::Symbol)
    value in LINE_PARAMETER_BASES || throw(
        ArgumentError("basis must be :per_length or :total; got :$value"),
    )
    return value
end

"""
    SeriesImpedance{T, Basis}

Store a square series-impedance matrix over frequency. Values use \\[Ω/m\\]
when `Basis` is `:per_length` and \\[Ω\\] when it is `:total`.
"""
struct SeriesImpedance{T, Basis} <: AbstractArray{T, 3}
    "Complex series-impedance tensor with dimensions conductor × conductor × frequency."
    values::Array{T, 3}

    function SeriesImpedance{T, Basis}(values::Array{T, 3}) where {T, Basis}
        Basis isa Symbol || throw(
            ArgumentError("basis must be :per_length or :total; got $(repr(Basis))"),
        )
        _check_basis(Basis)
        return new{T, Basis}(values)
    end
end

"""
    ShuntAdmittance{T, Basis}

Store a square shunt-admittance matrix over frequency. Values use \\[S/m\\]
when `Basis` is `:per_length` and \\[S\\] when it is `:total`.
"""
struct ShuntAdmittance{T, Basis} <: AbstractArray{T, 3}
    "Complex shunt-admittance tensor with dimensions conductor × conductor × frequency."
    values::Array{T, 3}

    function ShuntAdmittance{T, Basis}(values::Array{T, 3}) where {T, Basis}
        Basis isa Symbol || throw(
            ArgumentError("basis must be :per_length or :total; got $(repr(Basis))"),
        )
        _check_basis(Basis)
        return new{T, Basis}(values)
    end
end

function SeriesImpedance(A::AbstractArray{T, 3}; basis::Symbol = :per_length) where {T}
    _check_basis(basis)
    return SeriesImpedance{T, basis}(Array(A))
end

function ShuntAdmittance(A::AbstractArray{T, 3}; basis::Symbol = :per_length) where {T}
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

`Basis` is either `:per_length` or `:total`. Per-length values are stored in
Ω/m and S/m; total values are stored in Ω and S. Frequencies are stored in Hz.
"""
struct LineParameters{
    T <: COMPLEXSCALAR,
    U <: REALSCALAR,
    D <: LineParamsDomain,
    Basis
}
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
            T <: COMPLEXSCALAR,
            U <: REALSCALAR,
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
) where {T <: COMPLEXSCALAR, U <: REALSCALAR, Basis}
    return LineParameters(PhaseDomain, Z, Y, f)
end

function LineParameters(
        Z::SeriesImpedance{TZ, ZBasis},
        Y::ShuntAdmittance{TY, YBasis},
        f::AbstractVector{U}
) where {
        TZ <: COMPLEXSCALAR,
        TY <: COMPLEXSCALAR,
        U <: REALSCALAR,
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
        basis::Symbol = :per_length
) where {
        D <: LineParamsDomain,
        TZ <: COMPLEXSCALAR,
        TY <: COMPLEXSCALAR,
        U <: REALSCALAR
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
        basis::Symbol = :per_length
) where {
        TZ <: COMPLEXSCALAR,
        TY <: COMPLEXSCALAR,
        U <: REALSCALAR
}
    return LineParameters(PhaseDomain, Z, Y, f; basis)
end

@inline domain(::Type{<:LineParameters{T, U, D}}) where {T, U, D <: LineParamsDomain} = D
@inline domain(lp::LineParameters) = domain(typeof(lp))
@inline basis(::Type{<:LineParameters{T, U, D, Basis}}) where {T, U, D, Basis} = Basis
@inline basis(::LineParameters{T, U, D, Basis}) where {T, U, D, Basis} = Basis

frequencies(lp::LineParameters) = lp.f
nconductors(lp::LineParameters) = size(lp.Z, 1)
nfrequencies(lp::LineParameters) = length(lp.f)

series_impedance(lp::LineParameters) = lp.Z
shunt_admittance(lp::LineParameters) = lp.Y
series_impedance(impedance::SeriesImpedance) = impedance
shunt_admittance(admittance::ShuntAdmittance) = admittance
"""
    Z(parameters[, i, j[, k]])
    Y(parameters[, i, j[, k]])

Return series impedance or shunt admittance. With `(i, j)`, return the complete
frequency response at that matrix position. `k` may be one index, a range, or
`:`. Stored units follow [`basis`](@ref): \\[Ω/m\\] and \\[S/m\\] for
`:per_length`, or \\[Ω\\] and \\[S\\] for `:total`.
"""
Z(lp::LineParameters) = lp.Z
Y(lp::LineParameters) = lp.Y

Z(impedance::SeriesImpedance) = impedance.values
Z(impedance::SeriesImpedance, i, j) = view(impedance.values, i, j, :)
Z(impedance::SeriesImpedance, i, j, k) = impedance.values[i, j, k]
Y(admittance::ShuntAdmittance) = admittance.values
Y(admittance::ShuntAdmittance, i, j) = view(admittance.values, i, j, :)
Y(admittance::ShuntAdmittance, i, j, k) = admittance.values[i, j, k]

R(impedance::SeriesImpedance, args...) = real.(Z(impedance, args...))
X(impedance::SeriesImpedance, args...) = imag.(Z(impedance, args...))
G(admittance::ShuntAdmittance, args...) = real.(Y(admittance, args...))
B(admittance::ShuntAdmittance, args...) = imag.(Y(admittance, args...))
resistance(impedance::SeriesImpedance, args...) = R(impedance, args...)
reactance(impedance::SeriesImpedance, args...) = X(impedance, args...)
conductance(admittance::ShuntAdmittance, args...) = G(admittance, args...)
susceptance(admittance::ShuntAdmittance, args...) = B(admittance, args...)

@inline Z(lp::LineParameters, i, j) = view(lp.Z.values, i, j, :)
@inline Y(lp::LineParameters, i, j) = view(lp.Y.values, i, j, :)
@inline Z(lp::LineParameters, i, j, k) = lp.Z.values[i, j, k]
@inline Y(lp::LineParameters, i, j, k) = lp.Y.values[i, j, k]

R(lp::LineParameters, args...) = real.(Z(lp, args...))
X(lp::LineParameters, args...) = imag.(Z(lp, args...))
G(lp::LineParameters, args...) = real.(Y(lp, args...))
B(lp::LineParameters, args...) = imag.(Y(lp, args...))

@inline function _angular_frequencies(lp::LineParameters, k)
    selected = lp.f[k]
    any(iszero, selected isa Number ? (selected,) : selected) && throw(
        DomainError(selected, "L and C are undefined at zero frequency"),
    )
    return 2π .* selected
end

"""
$(TYPEDSIGNATURES)

Return series inductance using the same frequency-selection grammar as
[`Z`](@ref), evaluated as:

```math
L(f) = \\frac{\\operatorname{Im} Z(f)}{2\\pi f}.
```

Units are \\[H/m\\] for `:per_length` and \\[H\\] for `:total`.

# Errors

Throws `DomainError` when a selected frequency is zero.
"""
function L(lp::LineParameters)
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

Units are \\[F/m\\] for `:per_length` and \\[F\\] for `:total`.

# Errors

Throws `DomainError` when a selected frequency is zero.
"""
function C(lp::LineParameters)
    any(iszero, lp.f) && throw(DomainError(lp.f, "C is undefined at zero frequency"))
    return imag.(lp.Y.values) ./ reshape(2π .* lp.f, 1, 1, :)
end

L(lp::LineParameters, i, j) = L(lp, i, j, :)
C(lp::LineParameters, i, j) = C(lp, i, j, :)
L(lp::LineParameters, i, j, k) = imag.(Z(lp, i, j, k)) ./ _angular_frequencies(lp, k)
C(lp::LineParameters, i, j, k) = imag.(Y(lp, i, j, k)) ./ _angular_frequencies(lp, k)

resistance(lp::LineParameters, args...) = R(lp, args...)
reactance(lp::LineParameters, args...) = X(lp, args...)
inductance(lp::LineParameters, args...) = L(lp, args...)
conductance(lp::LineParameters, args...) = G(lp, args...)
susceptance(lp::LineParameters, args...) = B(lp, args...)
capacitance(lp::LineParameters, args...) = C(lp, args...)

_line_component_values(::Val{:R}, object, frequencies) = R(object)
_line_component_values(::Val{:X}, object, frequencies) = X(object)
_line_component_values(::Val{:G}, object, frequencies) = G(object)
_line_component_values(::Val{:B}, object, frequencies) = B(object)
_line_component_values(::Val{:Z_re}, object, frequencies) = R(object)
_line_component_values(::Val{:Z_im}, object, frequencies) = X(object)
_line_component_values(::Val{:Z_abs}, object, frequencies) = abs.(Z(object))
function _line_component_values(::Val{:Z_angle}, object, frequencies)
    angle.(Z(object)) .* (180 / π)
end
_line_component_values(::Val{:Y_re}, object, frequencies) = G(object)
_line_component_values(::Val{:Y_im}, object, frequencies) = B(object)
_line_component_values(::Val{:Y_abs}, object, frequencies) = abs.(Y(object))
function _line_component_values(::Val{:Y_angle}, object, frequencies)
    angle.(Y(object)) .* (180 / π)
end

_line_component_values(::Val{:L}, parameters::LineParameters, frequencies) = L(parameters)
_line_component_values(::Val{:C}, parameters::LineParameters, frequencies) = C(parameters)

function _frequency_component_values(component::Symbol, reactive, frequencies)
    any(iszero, frequencies) && throw(
        DomainError(frequencies, "$component is undefined at zero frequency"),
    )
    return reactive ./ reshape(2π .* frequencies, 1, 1, :)
end

function _line_component_values(::Val{:L}, object::SeriesImpedance, frequencies)
    _frequency_component_values(:L, X(object), frequencies)
end
function _line_component_values(::Val{:C}, object::ShuntAdmittance, frequencies)
    _frequency_component_values(:C, B(object), frequencies)
end
