defaults(::Type{Formula{:bottommost}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the deepest of the described `layers` by the formula `Formula{:bottommost}`. Its equation is
in [`EarthProperties`](@ref).

# Arguments

- `Formula{:bottommost}`: model of the equivalent earth.
- `layers`: the layers below the air, as the law describes them at `f`.
- `f`: frequency \\[Hz\\].

# Returns

- The last of `layers`.
"""
function (::Formula{:bottommost})(layers::Tuple{Vararg{EarthLayer{T}}}, f::T) where {T <:
                                                                                     Real}
    return last(layers)
end
