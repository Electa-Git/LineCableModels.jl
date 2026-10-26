"""
$(TYPEDSIGNATURES)

Return the distinct formulas of `formulation` in slot order. An earth formula that fills both
earth slots occurs once.
"""
function formulas(formulation::CoaxialFormulation{<:InternalImpedance.Formula,
        <:InsulationImpedance.Formula, <:InsulationAdmittance.Formula,
        <:SemiconAdmittance.Formula, E, E}) where {E <: EarthReturn.Formula}
    return (
        formulation.internal, formulation.insulation, formulation.insulation_admittance,
        formulation.semicon_admittance, formulation.earth_impedance, formulation.shunt,
        formulation.resistivity, formulation.properties)
end

# Earth slots of different types each contribute their formula.
function formulas(formulation::CoaxialFormulation)
    return (
        formulation.internal, formulation.insulation, formulation.insulation_admittance,
        formulation.semicon_admittance, formulation.earth_impedance,
        formulation.earth_admittance, formulation.shunt, formulation.resistivity,
        formulation.properties)
end

"""
$(TYPEDSIGNATURES)

Return the admittivity of `material` at `f` \\[Hz\\] in \\[S/m\\], by the law of its
material class in `formulation`.
"""
function admittivity(
        formulation::CoaxialFormulation, material::Material{
            T, :insulator}, f) where {T <: Real}
    return admittivity(formulation.insulation_admittance, material, f)
end

# A semicon layer takes the semicon admittivity law.
function admittivity(
        formulation::CoaxialFormulation, material::Material{
            T, :semicon}, f) where {T <: Real}
    return admittivity(formulation.semicon_admittance, material, f)
end

"""
$(TYPEDSIGNATURES)

Return `problem` unchanged for a formula without a condition on the problem.
"""
validate(problem::LineParametersProblem, ::AbstractFormulation) = problem

"""
$(TYPEDSIGNATURES)

Build the arrays for the coaxial computation of `problem` once, then add the arrays of each
distinct formula of `formulation`. Fill `gamma` with the longitudinal wavenumber \\[1/m\\]
at each frequency. The matrices store impedances \\[Ω/m\\], potential coefficients \\[m/F\\]
and their reductions. `rho` stores the resistivities of the conductors \\[Ω·m\\], and `kappa`
the admittivities of the dielectric layers \\[S/m\\].
"""
function initialize_buffers(
        formulation::CoaxialFormulation, problem::LineParametersProblem,
        input::CoaxialInput{T}) where {T <: Real}
    n = length(input.cable)
    l = length(input.layer_material)
    c = length(input.radius)
    m = length(input.kept)
    e = length(input.eliminated)
    nf = length(problem.frequencies)
    base = (rho = Vector{T}(undef, n),
        kappa = Vector{Complex{T}}(undef, l),
        gamma = Vector{Complex{T}}(undef, nf),
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
        formulas(formulation); init = base)
end

"""
$(TYPEDSIGNATURES)

Write the series impedance `Z` \\[Ω/m\\] and shunt admittance `Y` \\[S/m\\] at each
frequency of `problem`. Evaluate the conductor, dielectric and earth formulas, assemble the
impedance and potential matrices, and apply the selected reductions. The admittance is
``Y = s P^{-1}``, where ``P`` is the reduced potential matrix \\[m/F\\] and
``s = \\mathrm{j} 2 \\pi f`` \\[1/s\\]. When `options.trace`, copy the matrices inside
the cables and of the earth return into `trace` before assembly. Return `Z`.
"""
function solve!(
        Z::AbstractArray{Complex{T}, 3}, Y::AbstractArray{Complex{T}, 3}, trace,
        problem::LineParametersProblem{T}, input::CoaxialInput{T},
        formulation::CoaxialFormulation, buffers,
        options::ComputationOptions
) where {T <: Real}
    (; rho, kappa, Zin, Pin, Zg, Pg, reduced, factor, coupling, rhs, identity) = buffers
    for i in eachindex(rho)
        rho[i] = resistivity(formulation.resistivity, input.material[i], problem.temperature)
    end
    for (k, f) in enumerate(problem.frequencies)
        Zk = view(Z,:,:,k)
        Yk = view(Y,:,:,k)
        s = Complex{T}(im * (2 * π * f))
        for l in eachindex(kappa)
            kappa[l] = admittivity(formulation, input.layer_material[l], f)
        end
        media = map(layer -> properties(formulation.properties, layer, f), problem.earth.layers)
        impedance!(formulation.internal, formulation.insulation, Zin, input, rho, s)
        potential!(formulation.shunt, Pin, input, kappa, s)
        if options.trace
            copyto!(view(trace.Zin,:,:,k), Zin)
            copyto!(view(trace.Pin,:,:,k), Pin)
        end
        # The map compiles each call over the tuple of formulas, where foreach's fold dispatches at run time.
        map(formula -> prepare!(formula, buffers, input, media, s, k),
            formulas(formulation))
        impedance!(formulation.earth_impedance, Zg, buffers)
        potential!(formulation.earth_admittance, Pg, buffers)
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

Compute the line parameters of `problem` by the coaxial `formulation` on `LCMBackend`.
Validate the inputs, build the numerical input and buffers, and solve at each frequency.
Return the series impedances \\[Ω/m\\] and shunt admittances \\[S/m\\] in a
[`LineParameters`](@ref) record, with the retained conductor labels and the stage matrices
when `options.trace`.
"""
function compute(
        ::LCMBackend, problem::LineParametersProblem{T},
        formulation::CoaxialFormulation,
        options::ComputationOptions{LCMBackend}
) where {T <: Real}
    validate(problem, formulation)
    map(formula -> validate(problem, formula), formulas(formulation))
    input = CoaxialInput(problem, formulation.options)
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
