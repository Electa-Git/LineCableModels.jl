defaults(::Type{Formula{:series}}) = (;)

"""
$(TYPEDSIGNATURES)

Write into `Pin` the potential coefficients inside the cables, conductor by conductor, from the
series chain of their annular dielectric layers. A layer between the radii ``a`` and ``b``, of
admittivity ``\\kappa``, has the coefficient

```math
p = \\frac{s \\ln(b / a)}{2 \\pi \\kappa},
```

and zero when ``a`` is zero or equals ``b``. The coefficients of the layers in the interval of a
conductor add to its own coefficient. In each concentric assembly, the entry between two
conductors sums the own coefficients from the outer one of the two outward.

# Arguments

- `Pin`: matrix that receives the coefficients.
- `input`: numerical input of the cables.
- `kappa`: admittivity of each dielectric layer \\[S/m\\].
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].
- `coefficient`: buffer with one entry per dielectric layer.
- `tails`: buffer with one entry per conductor of the largest assembly.

# Returns

- `Pin`, its entries in \\[m/F\\].
"""
function potential!(::Formula{:series}, Pin::AbstractMatrix{Complex{T}},
        input::CoaxialInput{T}, kappa::AbstractVector{Complex{T}}, s::Complex{T},
        coefficient::AbstractVector{Complex{T}},
        tails::AbstractVector{Complex{T}}) where {T <: Real}
    fill!(Pin, zero(Complex{T}))
    for l in eachindex(coefficient)
        r_in = input.layer_r_in[l]
        r_ex = input.layer_r_ex[l]
        coefficient[l] = iszero(r_in) || r_in == r_ex ? zero(s) :
                         s / (2π * kappa[l] / log(r_ex / r_in))
    end
    for conductors in input.assemblies
        count = length(conductors)
        for position in count:-1:1
            i = conductors[position]
            own = zero(s)
            for l in input.layers[i]
                own += coefficient[l]
            end
            tails[position] = position == count ? own : own + tails[position + 1]
        end
        for row in 1:count, column in 1:count

            Pin[conductors[row], conductors[column]] += tails[max(row, column)]
        end
    end
    return Pin
end
