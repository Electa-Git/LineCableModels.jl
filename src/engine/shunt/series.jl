defaults(::Type{Formula{:series}}) = (;)

function (::Formula{:series})(
        r_in::T, r_ex::T, κ::Complex{T}, s::Complex{T}) where {T <: Real}
    return iszero(r_in) || r_in == r_ex ? zero(s) : s / (2π * κ / log(r_ex / r_in))
end

function (::Formula{:series})(
        n::Int, r_in::T, r_ex::T, κ::Complex{T}, s::Complex{T}) where {T <: Real}
    p = s / (2π * n * κ)
    ρ = (r_in / r_ex)^(2n)
    inner = p * (1 + ρ) / (1 - ρ)
    return (; inner, outer = inner, transfer = 2p * (r_in / r_ex)^n / (1 - ρ))
end

"""
$(TYPEDSIGNATURES)

Write into `Pin` the potential coefficients inside the cables, conductor by conductor, by
the formula `Formula{:series}`, whose equations are in the section "Formulas" of
[`Shunt`](@ref).

# Arguments

- `shunt`: annular dielectric series-chain formula.
- `Pin`: destination matrix of potential coefficients \\[m/F\\].
- `buffers`: computation arrays, with the image matrix of each interior.
- `input`: numerical input of the cables.
- `kappa`: admittivity of each dielectric layer \\[S/m\\].
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].

# Returns

- `Pin`, its entries in \\[m/F\\].
"""
function potential!(shunt::Formula{:series}, Pin::AbstractMatrix{Complex{T}}, buffers,
        input::CoaxialInput{T}, kappa::AbstractVector{Complex{T}},
        s::Complex{T}) where {T <: Real}
    fill!(Pin, zero(Complex{T}))
    for (a, conductors) in enumerate(input.assemblies)
        # The members' conductors precede the chain's from the first conductor of its first
        # member, and each charge inside a layer adds the layer's coefficient.
        start = isempty(input.members[a]) ? first(conductors) :
                first(first(buffers.interiors[a].ranges))
        count = length(conductors)
        running = zero(s)
        for position in count:-1:1
            i = conductors[position]
            own = zero(s)
            for l in input.dielectric_ranges[i]
                own += shunt(input.dielectric_r_in[l], input.dielectric_r_ex[l], kappa[l], s)
            end
            running = position == count ? own : own + running
            for row in start:i
                Pin[row, i] += running
            end
            for column in start:(i - 1)
                Pin[i, column] += running
            end
        end
        (isempty(input.members[a]) || isempty(conductors)) && continue
        # The gap between the members and the bore, and the returned fields, in the bore.
        f = first(input.interior_ranges[a])
        c = input.dielectric_r_ex[f]
        p_1 = shunt(1, zero(T), T(Inf), kappa[f], s).inner
        shell = running + sum(
            shunt(input.dielectric_r_in[l], input.dielectric_r_ex[l], kappa[l], s)
            for l in input.interior_ranges[a][2:end]; init = zero(s))
        (; potential, ranges, d) = buffers.interiors[a]
        for κ in eachindex(ranges), β in eachindex(ranges)

            term = shunt(d[β, κ], c, kappa[f], s) + shell + p_1 * potential[β, κ]
            for column in ranges[κ], row in ranges[β]

                Pin[row, column] += term
            end
        end
    end
    return Pin
end
