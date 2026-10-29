defaults(::Type{Formula{:series}}) = (;)

"""
$(TYPEDSIGNATURES)

Write into `Pin` the potential coefficients inside the cables, conductor by conductor, by
the formula `Formula{:series}`, whose equations are in the section "Formulas" of
[`Shunt`](@ref).

# Arguments

- `Formula{:series}`: annular dielectric series-chain formula.
- `Pin`: destination matrix of potential coefficients \\[m/F\\].
- `input`: numerical input of the cables.
- `kappa`: admittivity of each dielectric layer \\[S/m\\].
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].

# Returns

- `Pin`, its entries in \\[m/F\\].
"""
function potential!(::Formula{:series}, Pin::AbstractMatrix{Complex{T}},
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
                r_in = input.layer_r_in[l]
                r_ex = input.layer_r_ex[l]
                own += iszero(r_in) || r_in == r_ex ? zero(s) :
                       s / (2π * kappa[l] / log(r_ex / r_in))
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
