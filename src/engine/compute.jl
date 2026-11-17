"""
$(TYPEDSIGNATURES)

Return the admittivity of a medium from its resistivity and its relative permittivity at a
complex frequency:

```math
\\tilde{\\sigma} = \\frac{1}{\\rho} + s \\varepsilon_0 \\varepsilon_r,
```

where ``\\tilde{\\sigma}`` is the admittivity \\[S/m\\], ``\\rho`` the resistivity
\\[Ω·m\\], ``\\varepsilon_r`` the relative permittivity \\[dimensionless\\] and
``\\varepsilon_0`` the vacuum permittivity \\[F/m\\]. The complex frequency is
``s = \\mathrm{j} 2 \\pi f`` \\[1/s\\] at the frequency ``f`` \\[Hz\\].

# Arguments

- `rho`: resistivity \\[Ω·m\\].
- `eps_r`: relative permittivity \\[dimensionless\\].
- `s`: complex frequency \\[1/s\\].

# Returns

- Admittivity \\[S/m\\].
"""
function admittivity(rho::Real, eps_r::Real, s::Complex)
    return inv(rho) + s * (ε₀ * eps_r)
end

# An earth layer, declared or equivalent, takes the admittivity of its resistivity and its
# relative permittivity.
function admittivity(layer::Union{EarthLayer, EquivalentEarthLayer}, s::Complex)
    return admittivity(layer.rho, layer.eps_r, s)
end

# A conductor takes its resistivity at the temperature by the resistivity law of the
# formulation.
function admittivity(formulation::Union{CoaxialFormulation, FEM},
        material::Material{T, :conductor}, f, θ) where {T <: Real}
    return admittivity(formulation.resistivity(material, θ), material.eps_r,
        Complex{T}(im * (2 * π * f)))
end

"""
$(TYPEDSIGNATURES)

Evaluate the insulator admittivity law selected by the formulation, which reads the frequency
and not the temperature.

# Arguments

- `formulation`: coaxial or finite-element formulation with the material laws.
- `material`: insulator material.
- `f`: frequency \\[Hz\\].
- `θ`: temperature \\[°C\\].

# Returns

- Admittivity \\[S/m\\].
"""
function admittivity(formulation::Union{CoaxialFormulation, FEM},
        material::Material{T, :insulator}, f, θ) where {T <: Real}
    return formulation.insulation_admittance(material, f)
end

# A semicon layer takes the semicon admittivity law, which does not depend on the temperature.
function admittivity(formulation::Union{CoaxialFormulation, FEM},
        material::Material{T, :semicon}, f, θ) where {T <: Real}
    return formulation.semicon_admittance(material, f)
end

"""
$(TYPEDSIGNATURES)

Accept a problem for a formula without a condition on the problem.

# Arguments

- `problem`: line-parameters problem.
- `AbstractFormulation`: formulation or formula without a more specific validation method.

# Returns

- `problem`, unchanged.
"""
validate(problem::LineParametersProblem, ::AbstractFormulation) = problem

"""
$(TYPEDSIGNATURES)

Allocate the coaxial computation arrays, and fill the longitudinal wavenumbers and the
equivalent earth at each frequency. Add the arrays of each formula from the distinct formulas of
the formulation.

# Arguments

- `formulation`: formulas, longitudinal wavenumbers, and reductions.
- `problem`: line-parameters problem.
- `input`: numerical cable geometry and reduction indices.

# Returns

- A named tuple of conductor resistivities `rho` \\[Ω·m\\], dielectric admittivities `kappa`
  \\[S/m\\], and longitudinal wavenumbers `gamma` \\[1/m\\].
- The same tuple contains `earth`, the [`EquivalentEarthLayer`](@ref) of each frequency by the
  formulation's laws.
- The same tuple contains, for each medium, the complex conductivity `σ̃` \\[S/m\\] and the
  permeability `μ` \\[H/m\\].
- The same tuple contains `Zin`, `Zg` \\[Ω/m\\], `Pin`, `Pg` \\[m/F\\], reduction arrays, an
  identity matrix, and formula-specific arrays.
- The same tuple contains `interiors`, the record of each assembly's interior.

A record of `interiors` is empty where an assembly has no members. Otherwise, it holds the
arrays for the translations, sources, observations, reflections and returned fields.
The last axis of `reflections` holds the fields: 1 for the impedances, 2 for the potentials.
It also holds the factorization and its pivots, and the impedance and potential image
matrices. The key `d` holds the radius in \\[m\\] of each gap term of [`impedance!`](@ref).
The translations, sources and observations give the returned fields in closed form.

# Notes

A formula reads these terms and derives everything else itself:

| Term | Where a formula reads it | When the engine computes it |
| --- | --- | --- |
| `s` | the argument | the loop |
| geometry and layer map | `input`, with `x`, `y`, `radius`, `earth_layer`, the tubes and the dielectric layers | once, in `CoaxialInput` |
| `ρ_i` at the operating temperature, by the resistivity law | `buffers.rho[i]` | once, before the loop |
| `κ_l`, the admittivity of dielectric layer `l`, by its law | `buffers.kappa[l]` | per frequency, in `solve!` |
| `Γ_k` | `buffers.gamma[k]` | once, from the options |
| media as the law describes them and the model reduces them | `media[t]`, with `rho`, `eps_r`, `mu_r` | once, in `initialize_buffers` |
| `σ̃_t = 1/ρ_t + s ε₀ ε_r,t`, the [`admittivity`](@ref) of medium `t`, and `μ_t = μ₀ μ_r,t` | `buffers.σ̃[t]`, `buffers.μ[t]` | per frequency, by `to_buffer!` of the formulation |
"""
function initialize_buffers(
        formulation::CoaxialFormulation, problem::LineParametersProblem,
        input::CoaxialInput{T}) where {T <: Real}
    n = length(input.cable)
    l = length(input.dielectric_material)
    c = length(input.radius)
    nf = length(problem.frequencies)
    N = formulation.options.order
    # Each interior has its members and its bore as bodies. An assembly where no conductor
    # holds members has none.
    bodies = [isempty(input.members[a]) ? 0 : length(input.members[a]) + 1
              for a in eachindex(input.assemblies)]
    interiors = [(translations = zeros(T, 2N * B, 2N * B),
                     sources = zeros(T, 2N * B, max(B - 1, 0)),
                     observations = zeros(T, 2N * B, max(B - 1, 0)),
                     reflections = Array{Complex{T}, 4}(undef, 2, N, B, 2),
                     factor = Matrix{Complex{T}}(undef, 2N * B, 2N * B),
                     pivots = Vector{BlasInt}(undef, 2N * B),
                     returned = Matrix{Complex{T}}(undef, 2N * B, max(B - 1, 0)),
                     impedance = zeros(Complex{T}, max(B - 1, 0), max(B - 1, 0)),
                     potential = zeros(Complex{T}, max(B - 1, 0), max(B - 1, 0)),
                     d = zeros(T, max(B - 1, 0), max(B - 1, 0)))
                 for B in bodies]
    # The translations, sources and observations of the fields of the members and the bore,
    # in closed form. Each is a 2 × 2 block that acts on the cos and sin amplitudes of an order.
    # A binomial coefficient has the type of the nominal values of the radii.
    for (a, interior) in enumerate(interiors)
        isempty(input.members[a]) && continue
        (; translations, sources, observations, d) = interior
        P = bodies[a] - 1
        c_f = input.dielectric_r_ex[first(input.interior_ranges[a])]
        q = [complex(input.positions[m][1] - input.positions[a][1],
                 input.positions[m][2] - input.positions[a][2]) for m in input.members[a]]
        b = [input.r_ins_ex[last(input.assemblies[m])] for m in input.members[a]]
        # The radius of each gap term: a member's outer radius, or the distance between two axes.
        for (κ, n) in enumerate(input.members[a]), (β, m) in enumerate(input.members[a])

            d[β, κ] = m == n ? b[β] :
                      hypot(input.positions[m][1] - input.positions[n][1],
                input.positions[m][2] - input.positions[n][2])
        end
        for β in 1:bodies[a], ν in 1:N

            column = ((β - 1) * N + (ν - 1)) * 2 + 1
            for κ in 1:bodies[a], μ in 1:N

                row = ((κ - 1) * N + (μ - 1)) * 2 + 1
                if β == κ
                    continue
                elseif β <= P && κ <= P
                    δ = q[β] - q[κ]
                    t = binomial(oftype(nominal(c_f), ν + μ - 1), μ) * (b[β] / -δ)^ν *
                        (b[κ] / δ)^μ
                    translations[row:(row + 1), column:(column + 1)] .= [real(t) -imag(t);
                                                                         -imag(t) -real(t)]
                elseif β <= P && μ >= ν
                    t = binomial(oftype(nominal(c_f), μ - 1), ν - 1) * (b[β] / c_f)^ν *
                        (q[β] / c_f)^(μ - ν)
                    translations[row:(row + 1), column:(column + 1)] .= [real(t) -imag(t);
                                                                         imag(t) real(t)]
                elseif β > P && μ <= ν
                    t = binomial(oftype(nominal(c_f), ν), μ) * (q[κ] / c_f)^(ν - μ) *
                        (b[κ] / c_f)^μ
                    translations[row:(row + 1), column:(column + 1)] .= [real(t) imag(t);
                                                                         -imag(t) real(t)]
                end
            end
            for j in 1:P
                j == β && continue
                ζ = β <= P ? conj((b[β] / (q[j] - q[β]))^ν) : (q[j] / c_f)^ν
                observations[column:(column + 1), j] .= (real(ζ), imag(ζ))
                sources[column:(column + 1), j] .= (real(ζ) / ν, imag(ζ) / ν)
            end
        end
    end
    base = merge(
        (rho = Vector{T}(undef, n),
            kappa = Vector{Complex{T}}(undef, l),
            gamma = Vector{Complex{T}}(undef, nf),
            earth = formulation.equivalent_earth(problem, formulation.earth_properties),
            σ̃ = Vector{Complex{T}}(undef, 2),
            μ = Vector{T}(undef, 2),
            Zin = Matrix{Complex{T}}(undef, n, n),
            Pin = Matrix{Complex{T}}(undef, n, n),
            Zg = Matrix{Complex{T}}(undef, c, c),
            Pg = Matrix{Complex{T}}(undef, c, c),
            interiors),
        initialize_buffers(eliminate!, T, input.kept, input.eliminated))
    base.gamma .= formulation.options.Γ
    return foldl(
        (buffers, formula) -> initialize_buffers(formula, input, buffers),
        formulation.formulas; init = base)
end

"""
$(TYPEDSIGNATURES)

Return the impedance that a dielectric layer presents to a harmonic field at one radius
with a load at the other radius:

```math
Z' = z \\frac{Z + z \\tau}{z + Z \\tau}, \\qquad
\\tau = \\tanh\\left(n \\ln\\frac{r_\\mathrm{ex}}{r_\\mathrm{in}}\\right).
```

For the magnetic field, ``z = s \\mu_0 \\mu_r / (2\\pi n)`` \\[Ω/m\\] is the medium's
impedance from Schelkunoff 1934, equation (127) [Schelkunoff1934](@cite).
For the electric field, ``z = s / (2\\pi n \\kappa)`` \\[m/F\\] is its potential coefficient.
Here ``s`` is complex frequency \\[1/s\\], ``\\mu_r`` relative permeability \\[dimensionless\\],
and ``\\kappa`` admittivity \\[S/m\\]. The constant [`μ₀`](@ref) is vacuum permeability \\[H/m\\].

# Arguments

- `n`: positive harmonic order \\[dimensionless\\].
- `r_in`: positive inner radius \\[m\\].
- `r_ex`: positive outer radius \\[m\\].
- `z`: medium impedance \\[Ω/m\\] or potential coefficient \\[m/F\\] at order `n`.
- `Z`: load at the opposite radius, in the same units as `z`.

# Returns

- Impedance or potential coefficient at the layer's opposite radius, in the same
  units as `z`.

# Notes

The method applies a layer in either direction, inward or outward.
At zero width it returns the load, with outer-radius derivative ``n(z^2 - Z^2)/(zr)``
at the common radius ``r`` \\[m\\].
In the continuous limit of zero order, the magnetic result approaches
``Z + s \\mu_0 \\mu_r \\ln(r_\\mathrm{ex}/r_\\mathrm{in})/(2\\pi)``.
"""
function (::Union{AbstractInsulationImpedance, AbstractShunt})(
        n::Int, r_in::T, r_ex::T, z::Complex{T}, Z::Complex{T}) where {T <: Real}
    τ = tanh(n * log1p((r_ex - r_in) / r_in))
    return z * (Z + z * τ) / (z + Z * τ)
end

"""
$(TYPEDSIGNATURES)

Fill the terms that the engine provides to every formula at one frequency, one entry per medium
of `media`:

```math
\\tilde\\sigma_t = \\frac{1}{\\rho_t} + s \\varepsilon_0 \\varepsilon_{r,t}, \\qquad
\\mu_t = \\mu_0 \\mu_{r,t},
```

where ``\\rho_t`` is the resistivity \\[Ω·m\\], and ``\\varepsilon_{r,t}`` and ``\\mu_{r,t}`` are the
relative permittivity and permeability \\[dimensionless\\] of medium ``t``.

Above order 0, solve the image matrices of each interior at the orders 1 to `order`, for the
impedances and for the potentials, as [`impedance!`](@ref) gives them. The magnetic field uses
one outward fold per member, through each conductor's surface impedances and each dielectric
layer as a line. The bore folds them inward from the medium around the assembly. A conductor presents zero to the
electric field. For that field, each member folds its outer dielectric layers outward, and
the bore folds its interior layers inward. One solve per field gives its image matrix.

# Arguments

- `formulation`: coaxial formulation, with its formulas and the option `order`.
- `buffers`: computation arrays with `σ̃` \\[S/m\\], `μ` \\[H/m\\] and the interiors.
- `input`: numerical cable geometry.
- `media`: the air and the equivalent earth at the frequency.
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].
- `k`: frequency index.

# Returns

- `buffers`, with `σ̃`, `μ` and the image matrices of the interiors filled.
"""
function to_buffer!(
        formulation::CoaxialFormulation, buffers, input::CoaxialInput, media, s, k)
    # The air and the equivalent earth are of two types: the map compiles each medium's call.
    buffers.σ̃ .= map(medium -> admittivity(medium, s), media)
    buffers.μ .= map(medium -> μ₀ * medium.mu_r, media)
    N = formulation.options.order
    N == 0 && return buffers
    (; internal_impedance, insulation_impedance, shunt) = formulation
    (; rho, kappa) = buffers
    for (a, interior) in enumerate(buffers.interiors)
        isempty(input.members[a]) && continue
        bodies = length(input.members[a]) + 1
        (; translations, sources, observations,
            reflections, factor, pivots, returned) = interior
        U = 2N * bodies
        f = first(input.interior_ranges[a])
        μ_f = input.dielectric_material[f].mu_r
        # Beyond the assembly lies the filler of the interior around it, or free space.
        around = input.exterior_ranges[a]
        μ_around = isempty(around) ? oneunit(μ_f) :
                   input.dielectric_material[first(around)].mu_r
        # The magnetic field reflects each member outward and the bore inward.
        for n in 1:N
            z_n = s * μ₀ * μ_f / (2π * n)
            for (β, m) in enumerate(input.members[a])
                inner = input.interior_ranges[m]
                Z = s * μ₀ *
                    (isempty(inner) ? oneunit(μ_f) :
                     input.dielectric_material[first(inner)].mu_r) / (2π * n)
                for l in inner[2:end]
                    Z = insulation_impedance(
                        n, input.dielectric_r_in[l], input.dielectric_r_ex[l],
                        s * μ₀ * input.dielectric_material[l].mu_r / (2π * n), Z)
                end
                for i in input.assemblies[m]
                    surfaces = internal_impedance(n, input.r_in[i], input.r_ex[i], rho[i],
                        input.material[i].mu_r, s)
                    Z = surfaces.outer - surfaces.transfer^2 / (surfaces.inner + Z)
                    for l in input.dielectric_ranges[i]
                        Z = insulation_impedance(
                            n, input.dielectric_r_in[l], input.dielectric_r_ex[l],
                            s * μ₀ * input.dielectric_material[l].mu_r / (2π * n), Z)
                    end
                end
                reflections[1, n, β, 1] = reflections[2, n, β, 1] = (Z - z_n) / (Z + z_n)
            end
            Z = s * μ₀ * μ_around / (2π * n)
            for i in reverse(input.assemblies[a])
                for l in reverse(input.dielectric_ranges[i])
                    Z = insulation_impedance(
                        n, input.dielectric_r_in[l], input.dielectric_r_ex[l],
                        s * μ₀ * input.dielectric_material[l].mu_r / (2π * n), Z)
                end
                surfaces = internal_impedance(n, input.r_in[i], input.r_ex[i], rho[i],
                    input.material[i].mu_r, s)
                Z = surfaces.inner - surfaces.transfer^2 / (surfaces.outer + Z)
            end
            for l in reverse(input.interior_ranges[a][2:end])
                Z = insulation_impedance(
                    n, input.dielectric_r_in[l], input.dielectric_r_ex[l],
                    s * μ₀ * input.dielectric_material[l].mu_r / (2π * n), Z)
            end
            reflections[1, n, bodies, 1] = reflections[2, n, bodies, 1] = (Z - z_n) /
                                                                          (Z + z_n)
        end
        # The electric load is zero at a conductor. The fold applies its dielectric layers.
        for n in 1:N
            p_n = s / (2π * n * kappa[f])
            for (β, m) in enumerate(input.members[a])
                P = zero(s)
                for l in input.dielectric_ranges[last(input.assemblies[m])]
                    P = shunt(n, input.dielectric_r_in[l], input.dielectric_r_ex[l],
                        s / (2π * n * kappa[l]), P)
                end
                reflections[1, n, β, 2] = reflections[2, n, β, 2] = (P - p_n) / (P + p_n)
            end
            P = zero(s)
            for l in reverse(input.interior_ranges[a][2:end])
                P = shunt(n, input.dielectric_r_in[l], input.dielectric_r_ex[l],
                    s / (2π * n * kappa[l]), P)
            end
            reflections[1, n, bodies, 2] = reflections[2, n, bodies, 2] = (P - p_n) /
                                                                          (P + p_n)
        end
        for (reflection, image) in zip(eachslice(reflections; dims = 4),
            (interior.impedance, interior.potential))
            # The returned fields solve (I − R K) x = R t, one column per source member.
            for column in 1:U, row in 1:U

                factor[row, column] = (row == column) -
                                      reflection[row] * translations[row, column]
            end
            for j in axes(returned, 2), row in 1:U

                returned[row, j] = reflection[row] * sources[row, j]
            end
            ldiv!(factorize!(factor, pivots), returned)
            for j in axes(image, 2), i in axes(image, 1)

                image[i, j] = sum(observations[row, i] * returned[row, j] for row in 1:U)
            end
        end
    end
    return buffers
end

"""
$(TYPEDSIGNATURES)

Write the series impedances and shunt admittances at every frequency. Evaluate material laws
and cable formulas, and move the terms of each formula to its block. Assemble the conductor
matrices, and apply bundle and Kron reductions. The conductor matrices add the earth-return
matrices of the cables to the matrices inside them,

```math
Z_{ij} = Z_{\\mathrm{in}, ij} + Z_{g, c_i c_j}, \\qquad
P_{ij} = P_{\\mathrm{in}, ij} + P_{g, c_i c_j},
```

where ``c_i`` is the cable of conductor ``i``. The admittance is

```math
Y = s P_{\\mathrm{r}}^{-1}, \\qquad s = \\mathrm{j} 2 \\pi f,
```

where ``P_{\\mathrm{r}}`` is the reduced potential matrix \\[m/F\\] and ``f`` is frequency
\\[Hz\\]. Apply the ideal transposition after reduction when selected.

# Arguments

- `Z`: destination series-impedance array, indexed by row, column, and frequency \\[Ω/m\\].
- `Y`: destination shunt-admittance array in the same order \\[S/m\\].
- `trace`: destinations for `Zin`, `Zg` \\[Ω/m\\] and `Pin`, `Pg` \\[m/F\\] at each
  frequency.
- `problem`: cable system, earth, temperature, and frequencies.
- `input`: numerical cable geometry and reduction indices.
- `formulation`: formulas and formulation controls.
- `buffers`: preallocated material, matrix, and formula arrays.
- `options`: computation controls. The option `trace` selects copies of the intermediate
  matrices before assembly.

# Returns

- `Z`, with computed entries \\[Ω/m\\]. The method also fills `Y` and the selected trace
  arrays.
"""
function solve!(
        Z::AbstractArray{Complex{T}, 3}, Y::AbstractArray{Complex{T}, 3}, trace,
        problem::LineParametersProblem{T}, input::CoaxialInput{T},
        formulation::CoaxialFormulation, buffers,
        options::ComputationOptions
) where {T <: Real}
    (; rho, kappa, Zin, Pin, Zg, Pg, reduced, factor,
        coupling, rhs, identity, pivots) = buffers
    for i in eachindex(rho)
        rho[i] = formulation.resistivity(input.material[i], problem.temperature)
    end
    for (k, f) in enumerate(problem.frequencies)
        Zk = view(Z,:,:,k)
        Yk = view(Y,:,:,k)
        s = Complex{T}(im * (2 * π * f))
        for l in eachindex(kappa)
            kappa[l] = admittivity(
                formulation, input.dielectric_material[l], f, problem.temperature)
        end
        media = (first(problem.earth.layers), buffers.earth[k])
        to_buffer!(formulation, buffers, input, media, s, k)
        impedance!(formulation.internal_impedance,
            formulation.insulation_impedance, Zin, buffers, input, rho, s)
        potential!(formulation.shunt, Pin, buffers, input, kappa, s)
        if options.trace
            copyto!(view(trace.Zin,:,:,k), Zin)
            copyto!(view(trace.Pin,:,:,k), Pin)
        end
        # The map compiles each call over the tuple of formulas, where foreach's fold dispatches at run time.
        map(formula -> to_buffer!(formula, buffers, input, media, s, k),
            formulation.formulas)
        impedance!(formulation, Zg, buffers, input, media, s, k)
        potential!(formulation, Pg, buffers, input, media, s, k)
        if options.trace
            copyto!(view(trace.Zg,:,:,k), Zg)
            copyto!(view(trace.Pg,:,:,k), Pg)
        end
        for j in eachindex(input.cable), i in eachindex(input.cable)

            Zin[i, j] += Zg[input.cable[i], input.cable[j]]
            Pin[i, j] += Pg[input.cable[i], input.cable[j]]
        end
        bundle!(Zin, input.bundles)
        bundle!(Pin, input.bundles)
        eliminate!(Zk, Zin, input.kept, input.eliminated, factor, coupling, rhs, pivots)
        eliminate!(
            reduced, Pin, input.kept, input.eliminated, factor, coupling, rhs, pivots)
        ldiv!(Yk, factorize!(reduced, pivots), identity)
        Yk .*= s
        if formulation.options.transposition
            transposition!(Zk)
            transposition!(Yk)
        end
    end
    return Z
end

"""
$(TYPEDSIGNATURES)

Validate the problem for the formulation and its formulas, build numerical input and
buffers, and solve the coaxial line parameters at each frequency.

# Arguments

- `LCMBackend`: backend that executes the coaxial formulas.
- `problem`: cable system, earth, temperature, and frequencies.
- `formulation`: coaxial formulas and formulation controls.
- `options`: backend controls. The option `trace` retains intermediate matrices when `true`.

# Returns

- A [`LineParameters`](@ref) record with series impedances in \\[Ω/m\\], shunt admittances
  in \\[S/m\\], conductor labels, and the selected trace arrays.

# Errors

- Throws the errors of [`validate(problem, formulation)`](@ref validate),
  the [`CoaxialInput`](@ref) constructor, [`validate(input, formulation)`](@ref validate),
  and [`validate(input, formulation, problem)`](@ref validate).
"""
function compute(
        ::LCMBackend, problem::LineParametersProblem{T},
        formulation::CoaxialFormulation,
        options::ComputationOptions{LCMBackend}
) where {T <: Real}
    validate(problem, formulation)
    input = CoaxialInput(problem, formulation.options)
    validate(input, formulation)
    validate(input, formulation, problem)
    validate(merge(domain(problem, formulation.options), domain(input, problem)), formulation)
    buffers = initialize_buffers(formulation, problem, input)
    validate(buffers.earth, formulation.equivalent_earth, problem.frequencies)
    n = length(input.cable)
    c = length(input.radius)
    m = length(input.kept)
    nf = length(problem.frequencies)
    Z = Array{Complex{T}, 3}(undef, m, m, nf)
    Y = Array{Complex{T}, 3}(undef, m, m, nf)
    trace = options.trace ?
            (Zin = Array{Complex{T}, 3}(undef, n, n, nf),
        Pin = Array{Complex{T}, 3}(undef, n, n, nf),
        Zg = Array{Complex{T}, 3}(undef, c, c, nf),
        Pg = Array{Complex{T}, 3}(undef, c, c, nf)) :
            (Zin = Array{Complex{T}, 3}(undef, 0, 0, 0),
        Pin = Array{Complex{T}, 3}(undef, 0, 0, 0),
        Zg = Array{Complex{T}, 3}(undef, 0, 0, 0),
        Pg = Array{Complex{T}, 3}(undef, 0, 0, 0))
    solve!(Z, Y, trace, problem, input, formulation, buffers, options)
    return LineParameters{T}(Z, Y, problem.frequencies, input.terminals,
        convert(T, problem.system.line_length), formulation, trace)
end
