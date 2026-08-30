Base.eltype(::LineParameters{T}) where {T} = T
Base.eltype(::Type{LineParameters{T}}) where {T} = T

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
        selectdomain(lp.domain, selected),
        SeriesImpedance{T, Basis}(Array(view(lp.Z.values,:,:,selected))),
        ShuntAdmittance{T, Basis}(Array(view(lp.Y.values,:,:,selected))),
        selected_frequencies,
        lp.details
    )
end

"Return whether a scalar type carries explicit numerical uncertainty."
has_uncertainty_type(::Type) = false

function _result_unit(value, selector)
    quantity = Units.quantity(selector)
    return Units.native_unit(quantity, basis(value))
end
