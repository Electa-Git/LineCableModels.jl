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

# Arguments

- `internal`: formula of the surface impedances of each conductor.
- `insulation`: formula of the term of each dielectric interval.
- `Zin`: destination matrix of internal impedances \\[Ω/m\\].
- `input`: numerical input of the cables.
- `rho`: resistivity \\[Ω·m\\] of each conductor at the operating temperature.
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].

# Returns

- `Zin`, its entries in \\[Ω/m\\].
"""
function impedance!(internal::InternalImpedance.Formula,
        insulation::InsulationImpedance.Formula, Zin::AbstractMatrix{Complex{T}},
        input::CoaxialInput{T}, rho::AbstractVector{T}, s::Complex{T}) where {T <: Real}
    fill!(Zin, zero(Complex{T}))
    for conductors in input.assemblies
        inside = zero(s)
        for position in length(conductors):-1:1
            i = conductors[position]
            surfaces = internal(input.r_in[i], input.r_ex[i], rho[i],
                input.material[i].mu_r, s)
            term = insulation(input.r_ex[i], input.r_ins_ex[i],
                input.mu_r_ins[i], s)
            loop = surfaces.outer + inside + term
            for row in 1:(position - 1), column in 1:(position - 1)

                Zin[conductors[row], conductors[column]] += loop - 2 * surfaces.transfer
            end
            for row in 1:(position - 1)
                Zin[i, conductors[row]] += loop - surfaces.transfer
                Zin[conductors[row], i] += loop - surfaces.transfer
            end
            Zin[i, i] += loop
            inside = surfaces.inner
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
