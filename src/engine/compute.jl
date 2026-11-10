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
    base = merge(
        (rho = Vector{T}(undef, n),
            kappa = Vector{Complex{T}}(undef, l),
            gamma = Vector{Complex{T}}(undef, nf),
            earth = [formulation.equivalent_earth(
                         map(layer -> formulation.earth_properties(layer, f),
                             Base.tail(problem.earth.layers)),
                         f) for f in problem.frequencies],
            σ̃ = Vector{Complex{T}}(undef, 2),
            μ = Vector{T}(undef, 2),
            Zin = Matrix{Complex{T}}(undef, n, n),
            Pin = Matrix{Complex{T}}(undef, n, n),
            Zg = Matrix{Complex{T}}(undef, c, c),
            Pg = Matrix{Complex{T}}(undef, c, c)),
        initialize_buffers(eliminate!, T, input.kept, input.eliminated))
    base.gamma .= formulation.options.Γ
    return foldl(
        (buffers, formula) -> initialize_buffers(formula, input, buffers),
        formulation.formulas; init = base)
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

# Arguments

- `CoaxialFormulation`: coaxial formulation.
- `buffers`: computation arrays with `σ̃` \\[S/m\\] and `μ` \\[H/m\\].
- `input`: numerical cable geometry.
- `media`: the air and the equivalent earth at the frequency.
- `s`: complex frequency, `j 2π f`, in \\[1/s\\].
- `k`: frequency index.

# Returns

- `buffers`, with `σ̃` and `μ` filled.
"""
function to_buffer!(::CoaxialFormulation, buffers, input::CoaxialInput, media, s, k)
    # The air and the equivalent earth are of two types: the map compiles each medium's call.
    buffers.σ̃ .= map(medium -> admittivity(medium, s), media)
    buffers.μ .= map(medium -> μ₀ * medium.mu_r, media)
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
            formulation.insulation_impedance, Zin, input, rho, s)
        potential!(formulation.shunt, Pin, input, kappa, s)
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

- Throws `DimensionMismatch` when a longitudinal-wavenumber vector has a different length
  than the frequencies.
- Throws `DomainError` when a cable touches the interface for the earth-return formula
  `:martins2026`.
"""
function compute(
        ::LCMBackend, problem::LineParametersProblem{T},
        formulation::CoaxialFormulation,
        options::ComputationOptions{LCMBackend}
) where {T <: Real}
    validate(problem, formulation)
    input = CoaxialInput(problem, formulation.options)
    validate(input, formulation)
    validate(merge(domain(problem, formulation.options), domain(input)), formulation)
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
