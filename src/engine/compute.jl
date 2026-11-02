"""
$(TYPEDSIGNATURES)

Evaluate the insulator admittivity law selected by the formulation.

# Arguments

- `formulation`: coaxial or finite-element formulation with the material laws.
- `material`: insulator material.
- `f`: frequency \\[Hz\\].

# Returns

- Admittivity \\[S/m\\].
"""
function admittivity(
        formulation::Union{CoaxialFormulation, FEM}, material::Material{
            T, :insulator}, f) where {T <: Real}
    return formulation.insulation_admittance(material, f)
end

# A semicon layer takes the semicon admittivity law.
function admittivity(
        formulation::Union{CoaxialFormulation, FEM}, material::Material{
            T, :semicon}, f) where {T <: Real}
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

Allocate the coaxial computation arrays, fill the longitudinal wavenumbers, and add the
arrays of each formula from the distinct formulas of the formulation.

# Arguments

- `formulation`: formulas, longitudinal wavenumbers, and reductions.
- `problem`: line-parameters problem.
- `input`: numerical cable geometry and reduction indices.

# Returns

- A named tuple of conductor resistivities `rho` \\[Ω·m\\], dielectric admittivities `kappa`
  \\[S/m\\], and longitudinal wavenumbers `gamma` \\[1/m\\].
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
| media as the law describes them and the model reduces them | `media[t]`, with `rho`, `eps_r`, `mu_r`, `thickness` | per frequency, in `solve!` |
| `σ̃_t = 1/ρ_t + s ε₀ ε_r,t` and `μ_t = μ₀ μ_r,t` | `buffers.σ̃[t]`, `buffers.μ[t]` | per frequency, by `to_buffer!` of the formulation |
"""
function initialize_buffers(
        formulation::CoaxialFormulation, problem::LineParametersProblem,
        input::CoaxialInput{T}) where {T <: Real}
    n = length(input.cable)
    l = length(input.dielectric_material)
    c = length(input.radius)
    m = length(input.kept)
    e = length(input.eliminated)
    nf = length(problem.frequencies)
    base = (rho = Vector{T}(undef, n),
        kappa = Vector{Complex{T}}(undef, l),
        gamma = Vector{Complex{T}}(undef, nf),
        σ̃ = Vector{Complex{T}}(undef, 2),
        μ = Vector{T}(undef, 2),
        Zin = Matrix{Complex{T}}(undef, n, n),
        Pin = Matrix{Complex{T}}(undef, n, n),
        Zg = Matrix{Complex{T}}(undef, c, c),
        Pg = Matrix{Complex{T}}(undef, c, c),
        reduced = Matrix{Complex{T}}(undef, m, m),
        factor = Matrix{Complex{T}}(undef, e, e),
        coupling = Matrix{Complex{T}}(undef, m, e),
        rhs = Matrix{Complex{T}}(undef, e, m),
        identity = Matrix{Complex{T}}(I, m, m))
    base.gamma .= formulation.options.Γ
    return foldl(
        (buffers, formula) -> initialize_buffers(formula, problem, input, buffers),
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
    for t in eachindex(media)
        buffers.σ̃[t] = inv(media[t].rho) + s * (ε₀ * media[t].eps_r)
        buffers.μ[t] = μ₀ * media[t].mu_r
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
    (; rho, kappa, Zin, Pin, Zg, Pg, reduced, factor, coupling, rhs, identity) = buffers
    for i in eachindex(rho)
        rho[i] = formulation.resistivity(input.material[i], problem.temperature)
    end
    for (k, f) in enumerate(problem.frequencies)
        Zk = view(Z,:,:,k)
        Yk = view(Y,:,:,k)
        s = Complex{T}(im * (2 * π * f))
        for l in eachindex(kappa)
            kappa[l] = admittivity(formulation, input.dielectric_material[l], f)
        end
        air = first(problem.earth.layers)
        described = map(layer -> formulation.earth_properties(layer, f),
            Base.tail(problem.earth.layers))
        media = (air, formulation.equivalent_earth(described, f))
        to_buffer!(formulation, buffers, input, media, s, k)
        impedance!(formulation.internal_impedance,
            formulation.insulation_impedance, Zin, input, rho, s)
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
        eliminate!(Zk, Zin, input.kept, input.eliminated, factor, coupling, rhs)
        eliminate!(reduced, Pin, input.kept, input.eliminated, factor, coupling, rhs)
        ldiv!(Yk, lu!(reduced), identity)
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
    map(formula -> validate(problem, formula), formulation.formulas)
    input = CoaxialInput(problem, formulation.options)
    validate(input, formulation)
    buffers = initialize_buffers(formulation, problem, input)
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
