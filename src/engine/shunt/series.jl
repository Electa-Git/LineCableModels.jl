defaults(::Type{Formula{:series}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the potential coefficient of one annular dielectric layer by the formula
`Formula{:series}`. Its equations are in [`Shunt`](@ref).

# Arguments

- `Formula{:series}`: formula of the series chain.
- `r_in`: inner radius of the layer \\[m\\].
- `r_ex`: outer radius of the layer \\[m\\].
- `κ`: admittivity of the layer \\[S/m\\].
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].

# Returns

- Potential coefficient \\[m/F\\], zero when `r_in` is zero or equals `r_ex`.
"""
function (::Formula{:series})(
        r_in::T, r_ex::T, κ::Complex{T}, s::Complex{T}) where {T <: Real}
    return iszero(r_in) || r_in == r_ex ? zero(s) : s / (2π * κ / log(r_ex / r_in))
end

"""
$(TYPEDSIGNATURES)

Write into `Pin` the potential coefficients inside the cables, conductor by conductor, by
the formula `Formula{:series}`, whose equations are in the section "Formulas" of
[`Shunt`](@ref).

# Arguments

- `shunt`: annular dielectric series-chain formula.
- `Pin`: destination matrix of potential coefficients \\[m/F\\].
- `input`: numerical input of the cables.
- `kappa`: admittivity of each dielectric layer \\[S/m\\].
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].

# Returns

- `Pin`, its entries in \\[m/F\\].
"""
function potential!(shunt::Formula{:series}, Pin::AbstractMatrix{Complex{T}},
        input::CoaxialInput{T}, kappa::AbstractVector{Complex{T}},
        s::Complex{T}) where {T <: Real}
    fill!(Pin, zero(Complex{T}))
    for conductors in input.assemblies
        count = length(conductors)
        running = zero(s)
        for position in count:-1:1
            i = conductors[position]
            own = zero(s)
            for l in input.layers[i]
                own += shunt(input.layer_r_in[l], input.layer_r_ex[l], kappa[l], s)
            end
            running = position == count ? own : own + running
            for row in 1:position
                Pin[conductors[row], i] += running
            end
            for column in 1:(position - 1)
                Pin[i, conductors[column]] += running
            end
        end
    end
    return Pin
end
