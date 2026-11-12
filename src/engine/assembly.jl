"""
$(TYPEDSIGNATURES)

Write into `Zin` the series impedances inside the cables, conductor by conductor, by the
loop recursion of Ametani. In each concentric assembly, the current of a conductor returns
through the conductor outside it, from the outermost conductor inward. Conductor ``k`` has
the loop impedance

```math
z_k = z_{\\mathrm{outer}, k} + z_{\\mathrm{inner}, k + 1} + z_{\\mathrm{insulation}, k},
```

where ``z_{\\mathrm{outer}, k}`` is the outer surface impedance of conductor ``k`` and
``z_{\\mathrm{insulation}, k}`` is the term of its dielectric interval. The inner surface
impedance of the next conductor outward is ``z_{\\mathrm{inner}, k + 1}``, zero for the
outermost conductor. With the conductors of the assembly numbered ``1, \\ldots, n`` from the
inside out, the loop impedances and the transfer impedances ``z_{\\mathrm{transfer}, k}``
give the entries

```math
Z_{\\mathrm{in}, ij} = \\sum_{k = \\max(i, j)}^{n}
    \\left(z_k - \\left([i < k] + [j < k]\\right) z_{\\mathrm{transfer}, k}\\right),
```

where ``[\\cdot]`` is 1 when its condition holds and 0 otherwise. All impedance terms have
units \\[Ω/m\\].

An assembly whose bore holds members adds their loops, and each loop of its chain also includes
the currents of the members. Inside such a loop, the conductors begin at the first conductor of
the first member. Between the conductors of the members ``m`` and ``m'``, each with its members'
conductors, the entries add

```math
z_{mm'} = z_\\mathrm{insulation}(d_{mm'}, c) + \\sum_l z_{\\mathrm{insulation}, l}
    + z_{\\mathrm{inner}, \\beta} + \\frac{s \\mu_0 \\mu_f}{2 \\pi} A_{mm'},
```

where ``c`` is the outer radius of the filler, the region around the members, and ``\\mu_f``
its relative permeability. The radius ``d_{mm'}`` is the distance between the axes of two
members, and ``d_{mm}`` is the outer radius of the interval of ``m``. The sum runs over the
annuli between the filler and the bore ``\\beta``, whose inner surface impedance at order 0 is
``z_{\\mathrm{inner}, \\beta}``.

The image matrix ``A`` holds the fields that the members and the bore return at the orders 1
to `order` [Hoidalen2025](@cite). Each body reflects the field of the order ``n`` by
``R_n = (Z_n - z_n) / (Z_n + z_n)``. The impedance ``Z_n`` is the one it presents through its
own conductors and layers, and ``z_n`` is that of the filler. The translations ``\\mathbf K``
give each returned field on the other bodies, and the sources ``\\mathbf t`` give the field of
the members' currents. The values ``\\mathbf o`` give the returned fields at the members' axes,
and the returned fields ``\\mathbf x`` solve

```math
(\\mathbf I - \\mathbf R \\mathbf K) \\mathbf x = \\mathbf R \\mathbf t, \\qquad
\\mathbf A = \\mathbf o^\\mathsf{T} \\mathbf x,
```

where ``\\mathbf R`` is the diagonal of the reflections. At `order = 0` the image matrix is
zero.

# Arguments

- `internal`: formula of the surface impedances of each conductor.
- `insulation`: formula of the term of each dielectric interval.
- `Zin`: destination matrix of internal impedances \\[Ω/m\\].
- `buffers`: computation arrays, with the image matrix of each interior.
- `input`: numerical input of the cables.
- `rho`: resistivity \\[Ω·m\\] of each conductor at the operating temperature.
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].

# Returns

- `Zin`, its entries in \\[Ω/m\\].

# Notes

Above order 0, the medium beyond a cable is free space. A member with members of its own
presents its conductors and layers, with its filler inside them.
"""
function impedance!(internal::InternalImpedance.Formula,
        insulation::InsulationImpedance.Formula, Zin::AbstractMatrix{Complex{T}}, buffers,
        input::CoaxialInput{T}, rho::AbstractVector{T}, s::Complex{T}) where {T <: Real}
    fill!(Zin, zero(Complex{T}))
    for (a, conductors) in enumerate(input.assemblies)
        # The loops of the chain carry the currents of its members, whose conductors precede
        # the chain's from the first conductor of its first member.
        start = isempty(input.members[a]) ? first(conductors) :
                first(first(buffers.interiors[a].ranges))
        inside = zero(s)
        for position in length(conductors):-1:1
            i = conductors[position]
            surfaces = internal(input.r_in[i], input.r_ex[i], rho[i],
                input.material[i].mu_r, s)
            term = insulation(input.r_ex[i], input.r_ins_ex[i],
                input.mu_r_ins[i], s)
            loop = surfaces.outer + inside + term
            for row in start:(i - 1), column in start:(i - 1)

                Zin[row, column] += loop - 2 * surfaces.transfer
            end
            for row in start:(i - 1)
                Zin[i, row] += loop - surfaces.transfer
                Zin[row, i] += loop - surfaces.transfer
            end
            Zin[i, i] += loop
            inside = surfaces.inner
        end
        (isempty(input.members[a]) || isempty(conductors)) && continue
        # The members' outermost loops and their returned fields, in the bore of the chain.
        f = first(input.interior_ranges[a])
        c = input.dielectric_r_ex[f]
        μ_f = input.dielectric_material[f].mu_r
        z_1 = insulation(1, zero(T), T(Inf), μ_f, s).inner
        shell = inside + sum(
            insulation(input.dielectric_r_in[l], input.dielectric_r_ex[l],
                input.dielectric_material[l].mu_r, s)
            for l in input.interior_ranges[a][2:end]; init = zero(s))
        (; impedance, ranges, d) = buffers.interiors[a]
        for κ in eachindex(ranges), β in eachindex(ranges)

            term = insulation(d[β, κ], c, μ_f, s) + shell + z_1 * impedance[β, κ]
            for column in ranges[κ], row in ranges[β]

                Zin[row, column] += term
            end
        end
    end
    return Zin
end

"""
$(TYPEDSIGNATURES)

Write into `Zg` the earth-return impedances of the cables, pair by pair, from the formula of
each pair's domain in the slot `earth_impedance`. The recipe alone writes the matrix.

# Arguments

- `formulation`: coaxial formulation with its earth-return formulas.
- `Zg`: destination matrix of cable impedances \\[Ω/m\\].
- `buffers`: preallocated computation arrays.
- `input`: numerical cable geometry and connections.
- `media`: the air and the equivalent earth at the frequency.
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].
- `k`: frequency index.

# Returns

- `Zg`, its entries in \\[Ω/m\\].
"""
function impedance!(formulation::CoaxialFormulation, Zg::AbstractMatrix, buffers,
        input::CoaxialInput, media, s, k)
    return assemble!(Zg, formulation.earth_impedance, Val(:Z), buffers, input, media, s, k)
end

"""
$(TYPEDSIGNATURES)

Write into `Pg` the earth-return potential coefficients of the cables, pair by pair, from
the formula of each pair's domain in the slot `earth_admittance`. The recipe alone writes
the matrix.

# Arguments

- `formulation`: coaxial formulation with its earth-return formulas.
- `Pg`: destination matrix of cable potential coefficients \\[m/F\\].
- `buffers`: preallocated computation arrays.
- `input`: numerical cable geometry and connections.
- `media`: the air and the equivalent earth at the frequency.
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].
- `k`: frequency index.

# Returns

- `Pg`, its entries in \\[m/F\\].
"""
function potential!(formulation::CoaxialFormulation, Pg::AbstractMatrix, buffers,
        input::CoaxialInput, media, s, k)
    return assemble!(Pg, formulation.earth_admittance, Val(:P), buffers, input, media, s, k)
end

"""
$(TYPEDSIGNATURES)

Write into `M` the earth-return terms of the cables for the quantity of the tag `quantity`,
pair by pair. Each pair takes the formula of its domain in `slot`. A pair above the earth
takes the air formula, a pair in it the earth formula, and any other pair the mixed
formula. A diagonal entry takes the formula of its cable's domain.

# Arguments

- `M`: destination matrix of the quantity per unit length.
- `slot`: the air, earth and mixed formulas of one earth slot of the formulation.
- `quantity`: tag of the quantity, `Val(:Z)` for impedances \\[Ω/m\\] or `Val(:P)` for
  potential coefficients \\[m/F\\].
- `buffers`: preallocated computation arrays.
- `input`: numerical cable geometry and connections.
- `media`: the air and the equivalent earth at the frequency.
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].
- `k`: frequency index.

# Returns

- `M`, in the unit of the quantity.
"""
function assemble!(M::AbstractMatrix, slot, quantity::Val, buffers, input::CoaxialInput,
        media, s, k)
    (; air, earth, mixed) = slot
    for i in eachindex(input.radius)
        M[i, i] = (input.earth_layer[i] == 1 ? air : earth)(
            quantity, Val(:self), i, buffers, input, media, s, k)
    end
    for j in eachindex(input.radius), i in eachindex(input.radius)

        i == j && continue
        above_i, above_j = input.earth_layer[i] == 1, input.earth_layer[j] == 1
        M[i, j] = if above_i && above_j
            air(quantity, Val(:mutual), i, j, buffers, input, media, s, k)
        elseif !above_i && !above_j
            earth(quantity, Val(:mutual), i, j, buffers, input, media, s, k)
        else
            mixed(quantity, Val(:mutual), i, j, buffers, input, media, s, k)
        end
    end
    return M
end
