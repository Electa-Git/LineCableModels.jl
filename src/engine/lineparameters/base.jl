Base.eltype(::LineParameters{T}) where {T} = T
Base.eltype(::Type{LineParameters{T}}) where {T} = T

Base.size(value::SeriesImpedance) = size(value.values)
Base.size(value::SeriesImpedance, dimension::Int) = size(value.values, dimension)
Base.axes(value::SeriesImpedance) = axes(value.values)
Base.ndims(::Type{<:SeriesImpedance}) = 3
Base.eltype(::Type{SeriesImpedance{T}}) where {T} = T
Base.getindex(value::SeriesImpedance, indices...) = getindex(value.values, indices...)
Base.IndexStyle(::Type{<:SeriesImpedance}) = IndexCartesian()

Base.size(value::ShuntAdmittance) = size(value.values)
Base.size(value::ShuntAdmittance, dimension::Int) = size(value.values, dimension)
Base.axes(value::ShuntAdmittance) = axes(value.values)
Base.ndims(::Type{<:ShuntAdmittance}) = 3
Base.eltype(::Type{ShuntAdmittance{T}}) where {T} = T
Base.getindex(value::ShuntAdmittance, indices...) = getindex(value.values, indices...)
Base.IndexStyle(::Type{<:ShuntAdmittance}) = IndexCartesian()

function Base.getindex(
        lp::LineParameters{T},
        selector::Union{
            Integer, AbstractRange{<:Integer}, AbstractVector{<:Integer}, Colon}
) where {T}
    selected = selector isa Integer ? (selector:selector) : selector
    checkbounds(lp.f, selected)
    selected_frequencies = selector isa Integer ? lp.f[selector:selector] : lp.f[selected]
    return LineParameters(
        selectdomain(lp.domain, selected),
        SeriesImpedance{T}(Array(view(lp.Z.values,:,:,selected))),
        ShuntAdmittance{T}(Array(view(lp.Y.values,:,:,selected))),
        selected_frequencies,
        lp.terminals,
        lp.line_length,
        selectdetails(lp.details,lp.domain,selected)
    )
end

"""
Return whether a scalar type encodes explicit numerical uncertainty.
"""
has_uncertainty_type(::Type) = false
