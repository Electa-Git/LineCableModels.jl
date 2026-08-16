Base.eltype(::LineParametersProblem{T}) where {T} = T
Base.eltype(::Type{LineParametersProblem{T}}) where {T} = T
Base.eltype(::LineParameters{T}) where {T} = T
Base.eltype(::Type{LineParameters{T}}) where {T} = T
Base.eltype(::EMTWorkspace{T}) where {T} = T
Base.eltype(::Type{EMTWorkspace{T}}) where {T} = T

Base.size(value::SeriesImpedance) = size(value.values)
Base.size(value::SeriesImpedance, dimension::Int) = size(value.values, dimension)
Base.axes(value::SeriesImpedance) = axes(value.values)
Base.ndims(::Type{<:SeriesImpedance}) = 3
Base.eltype(::Type{SeriesImpedance{T, Basis}}) where {T, Basis} = T
Base.getindex(value::SeriesImpedance, indices...) = getindex(value.values, indices...)
Base.IndexStyle(::Type{<:SeriesImpedance}) = IndexCartesian()

Base.size(value::ShuntAdmittance) = size(value.values)
Base.size(value::ShuntAdmittance, dimension::Int) = size(value.values, dimension)
Base.axes(value::ShuntAdmittance) = axes(value.values)
Base.ndims(::Type{<:ShuntAdmittance}) = 3
Base.eltype(::Type{ShuntAdmittance{T, Basis}}) where {T, Basis} = T
Base.getindex(value::ShuntAdmittance, indices...) = getindex(value.values, indices...)
Base.IndexStyle(::Type{<:ShuntAdmittance}) = IndexCartesian()

function Base.getindex(
        lp::LineParameters{T, U, D, Basis},
        selector::Union{
            Integer, AbstractRange{<:Integer}, AbstractVector{<:Integer}, Colon}
) where {T, U, D, Basis}
    selected = selector isa Integer ? (selector:selector) : selector
    checkbounds(lp.f, selected)
    selected_frequencies = selector isa Integer ? lp.f[selector:selector] : lp.f[selected]
    return LineParameters(
        D,
        SeriesImpedance{T, Basis}(Array(view(lp.Z.values, :, :, selected))),
        ShuntAdmittance{T, Basis}(Array(view(lp.Y.values, :, :, selected))),
        selected_frequencies
    )
end

_has_uncertainty_type(::Type) = false

@inline function _basis_units(value, per_length_unit, total_unit)
    return basis(value) === :per_length ? per_length_unit : total_unit
end

function Base.show(io::IO, value::SeriesImpedance)
    print(
        io,
        "SeriesImpedance ",
        join(size(value), '×'),
        " [",
        _basis_units(value, "Ω/m", "Ω"),
        "]"
    )
end

function Base.show(io::IO, value::ShuntAdmittance)
    print(
        io,
        "ShuntAdmittance ",
        join(size(value), '×'),
        " [",
        _basis_units(value, "S/m", "S"),
        "]"
    )
end

function Base.show(io::IO, lp::LineParameters)
    element_type = eltype(lp)
    print(
        io,
        "LineParameters{$element_type} ",
        join(size(lp.Z), '×'),
        " [",
        basis(lp),
        "]"
    )
    _has_uncertainty_type(element_type) && print(io, " (±)")
end

function Base.show(io::IO, ::MIME"text/plain", lp::LineParameters)
    impedance_unit = _basis_units(lp, "Ω/m", "Ω")
    admittance_unit = _basis_units(lp, "S/m", "S")
    println(io, join(size(lp.Z), '×'), " LineParameters [", basis(lp), "]")
    println(io, "domain: ", nameof(domain(lp)), ", frequencies: ", nfrequencies(lp))
    println(io, "Z [$impedance_unit], first frequency slice:")
    show(io, MIME"text/plain"(), view(lp.Z.values, :, :, 1))
    print(io, "\nY [$admittance_unit], first frequency slice:\n")
    show(io, MIME"text/plain"(), view(lp.Y.values, :, :, 1))
end
