"""
$(TYPEDSIGNATURES)

Validate Unified's prescribed longitudinal wavenumber Γ \\[1/m\\]. A scalar
applies at every frequency. A nonempty vector follows the frequency order. Return `Γ`.
"""
function validate(Γ, ::Type{<:Union{EarthImpedance.Formula{:unified}, Formula{:unified}}})
    Γ isa Union{Number, AbstractVector} || throw(ArgumentError(
        "unified Γ must be a scalar or frequency-aligned vector [1/m]"))
    values = Γ isa Number ? (Γ,) : Γ
    !isempty(values) &&
    all(value -> value isa Number && !(value isa Bool) && isfinite(value), values) ||
        throw(ArgumentError("unified Γ must be a finite scalar or nonempty finite vector [1/m]"))
    return Γ
end

function formulation_options(
        owner::Type{<:Union{EarthImpedance.Formula{:unified}, Formula{:unified}}},
        options::FormulationOptions)
    argument = validate(get(options.data, :Γ, 0), owner)
    return FormulationOptions(merge(options.data,
        (; Γ = argument isa AbstractVector ? copy(argument) : argument)))
end

# The formula's construction validated Γ. Its expressions take the value as supplied.
function formulation_options(
        expression::Expression{<:Union{
            EarthImpedance.Formula{:unified}, Formula{:unified}}},
        ::Val{:Γ}, default, supplied)
    return supplied
end

function description(::Type{<:Formula{:unified}}; compact::Bool = false)
    compact ? "Unified" :
    "Unified circumferential earth potential with complete enclosed-current normalization"
end

function description(::Type{<:Formula{:unified}}, ::Val{:Γ}, value::Number; compact::Bool = false)
    "Γ="*string(value)*" m⁻¹"
end
function description(::Type{<:Formula{:unified}}, ::Val{:Γ}, value::AbstractVector; compact::Bool = false)
    "Γ=["*join(value, ", ")*"] m⁻¹ (frequency order)"
end

function computation_type(::Type{T},
        selected::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        frequencies) where {T <: Real}
    argument = selected.options.data.Γ
    argument isa AbstractVector && length(argument) != length(frequencies) &&
        throw(DimensionMismatch("unified Γ must contain one value per frequency sample"))
    argument isa Number &&
        return promote_type(T, typeof(real(argument)), typeof(imag(argument)))
    return foldl(argument; init = T) do scalar, value
        promote_type(scalar, typeof(real(value)), typeof(imag(value)))
    end
end

"""
$(TYPEDSIGNATURES)

Calculate both source coefficients of a conductor pair for circumferentially averaged fields
in two homogeneous half-spaces: the axial-field coefficient per source current \\[Ω/m\\] and
the source-charge coefficient of conductor voltage \\[m/F\\]. It takes the Unified formula of
either earth family. The exp(jωt) convention and caller-prescribed Γ \\[1/m\\]
give

```math
\\widetilde K_{ij}=\\mathcal Z_{ij}-\\Gamma^2\\mathcal P_{\\phi,ij}/(j\\omega).
```

The electric scalar-potential contribution is retained when Γ is nonzero. Air targets use
the interface z=0 as voltage reference. Buried targets use deep earth. With source
amplitudes q̃ and conductor voltages U,

```math
U=\\widetilde P\\widetilde q,\\qquad P_e T_I=\\widetilde P.
```

The source-potential coefficient includes the complete-field voltage relation, not just
scalar potential. The air endpoint and conductor field are combined before quadrature to
preserve cancellation. Source columns include exp(abs(real(κⱼrⱼ))) scaling, shared by both
coefficients and the enclosed-current matrices, and it cancels in the physical matrix solve.
`functor.input.pair` retains source-target geometry \\[m\\]. `functor.state` contains the
evaluated media of the system, and the workspace buffers store the circumferential factors.
The complete current map converts the axial-field coefficients to physical series impedance.
The direct-image term is evaluated once for both coefficients. Quadrature estimates remain
diagnostic warnings.

# Returns

- `(axial, potential)`: the scaled axial-field coefficient \\[Ω/m\\] and the scaled
  source-potential coefficient \\[m/F\\].

# Reference

User-supplied manuscript, *Unified circumferentially averaged framework for overhead,
buried, and mixed conductor systems*: complete-field current relation, voltage and
source-charge maps.
"""
function source_coefficients(
        formula::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        kind::Union{Val{:self}, Val{:mutual}}, source::Union{Val{1}, Val{2}},
        target::Union{Val{1}, Val{2}}, functor, workspace)
    u=functor.state
    pair=functor.input.pair
    buffers=workspace.buffers
    row, column=pair.row, pair.column
    hp, hq=abs(pair.heights[2]), abs(pair.heights[1])
    y=pair.separation
    r=workspace.plan.geometry.radius[row]
    average=buffers.circumference_average[row]
    argument=buffers.radial_argument[row]
    sp, sq=buffers.source_logscale[row], buffers.source_logscale[column]
    πT=one(u.jω)*π
    integration=functor.input.options.data.integration
    context=(
        formula = :unified, frequency = imag(u.jω)/(2π), receiver = row, source = column)
    direct=earth_direct(formula, kind, source, target, u, pair, r, average, argument, sp, sq)
    # The axial field, before the source potential: the integral records keep this order.
    medium=target === Val(1) ? 1 : 2
    z=u.jω/πT*average*earth_spectral_term(formula, Val(:Z), target, source, u,
        hp, hq, y, zero(r), sp+sq,
        integration.method, integration.options, buffers; context)
    z+=u.jω*u.mu[medium]/(2πT)*direct
    phi=zero(z)
    if !iszero(u.Γ)
        phi=u.jω/πT*average*earth_spectral_term(formula, Val(:phi), target, source, u,
            hp, hq, y, zero(r), sp+sq,
            integration.method, integration.options, buffers; context)
        phi+=u.jω/(2πT*u.sh[medium])*direct
    end
    axial=z-u.Γ^2/u.jω*phi
    # A buried target takes its voltage from deep earth.
    if target === Val(2)
        value=u.jω/(one(u.jω)*π)*average*earth_spectral_term(formula, Val(:voltage), Val(2),
            source, u, hp, hq, y, zero(r), sp+sq, integration.method,
            integration.options, buffers; context)
        return (axial, value+u.jω/(2*(one(u.jω)*π)*u.sh[2])*direct)
    end
    # An air target takes its voltage from the interface.
    if abs(real(nominal(argument)))<300
        R=typeof(float(nominal(real(u.jω))))
        padding=hq/2
        angle=min(R(π)/6, atan(R(nominal(hq))/(4max(R(nominal(y+r)), eps(R)))))
        angle=earth_contour_angle(u, angle)
        points=earth_spectral_points!(buffers.earth_spectrum, u, hq-padding, y, r, angle)
        g=(hp, hq, radius = r, padding, logscale = sq, i0minus = bessel_i0m1(argument))
        S=source === Val(1) ? 1 : 2
        kernel=AirVoltageSpectrum{S, typeof(u), typeof(g)}(u, g)
        scale=max(R(abs(nominal(u.k[2]))), inv(R(nominal(hq))))
        contour=scale*cis(angle)
        integral=SpectralIntegral(t->contour*earth_spectral_value(
            kernel, contour*t, hq-padding, y, zero(r)))
        points ./= scale
        push!(points, one(scale))
        value,
        _=integrate(integral, integration.method, integration.options,
            buffers; points, coordinate_type = R,
            context = merge(context, (term = :air_voltage,)), observations = buffers.observations)
        return (axial, u.jω/πT*value+u.jω/(2πT*u.sh[1])*direct)
    end
    value=u.jω/πT*average*earth_spectral_term(formula, Val(:voltage), Val(1), source, u,
        hp, hq, y, zero(r), sp+sq, integration.method, integration.options, buffers; context)
    value+=u.jω/(2πT*u.sh[1])*direct
    return (axial, value+u.jω/πT*earth_spectral_term(formula, Val(:air_reference), Val(1),
        source, u, zero(r), hq, y, r, sq, integration.method, integration.options, buffers;
        context))
end

function Expression(formula::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        pair::EarthPair)
    return Expression(formula, source_coefficients,
        Val(pair.row == pair.column ? :self : :mutual), Val.(layer_index(pair))...)
end

function formulation_options(::Expression{
        <:Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        typeof(source_coefficients),
        A}) where {
        A <: Tuple{
        Union{Val{:self}, Val{:mutual}}, Union{Val{1}, Val{2}}, Union{Val{1}, Val{2}}}}
    return FormulationOptions((Γ = 0, integration = (method = :quad, options = (;))))
end

function validate(reduction::EquivalentHomogeneous.Formula{:bottommost},
        ::Expression{<:Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
            typeof(source_coefficients)})
    return reduction
end

# The unified calculation computes the whole system, every conductor pair, even when its slot
# publishes some of them, and both coefficients of each pair: either output needs both. Its
# formulation options are common to all pairs, and the exterior circumferences must lie in one
# half-space and must not overlap.
function EarthPlan(formula::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        earth, model::EarthModel, physical::AbstractVector{<:EarthPair}, indices,
        geometry::NamedTuple)
    whole=invoke(EarthPlan,
        Tuple{Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
            Any, EarthModel, AbstractVector{<:EarthPair}, Any, NamedTuple},
        formula, earth, model, physical, eachindex(physical), geometry)
    calculation=only(whole.calculations)
    options=first(calculation.parts).options
    all(part->isequal(part.options, options), calculation.parts) ||
        throw(ArgumentError("the unified earth-return calculation requires common formulation options for all conductor pairs"))
    for entry in calculation.pairs
        pair=entry.pair
        target_radius=geometry.radius[pair.row]
        source_radius=geometry.radius[pair.column]
        # Physical shapes can be disjoint while their equivalent circles overlap.
        # The field coefficients integrate these circles, so their applicability
        # is checked here using the engine's geometry.
        if pair.row==pair.column
            target_radius<abs(pair.heights[1]) || throw(DomainError(target_radius,
                "each exterior circumference must lie wholly in one half-space"))
        else
            hypot(pair.separation, pair.heights[1]-pair.heights[2])>
            target_radius+source_radius || throw(DomainError((pair.row, pair.column),
                "exterior circumferences must not overlap"))
        end
    end
    published=(; formula, pairs = collect(indices))
    impedance=formula isa EarthImpedanceFormulation ? published : nothing
    admittance=formula isa EarthAdmittanceFormulation ? published : nothing
    return EarthPlan((merge(calculation, (; impedance, admittance)),))
end

# Unified's arithmetic reads the pair's geometry and the conductor radii, not its indices.
function same_physical_state(::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        a::EarthPair, b::EarthPair, geometry::NamedTuple)
    inputs(pair)=(pair.row==pair.column, pair.layers, pair.heights, pair.separation,
        geometry.radius[pair.row], geometry.radius[pair.column])
    return same_physical_state(inputs(a), inputs(b))
end

# An impedance and an admittance formula of Unified share one calculation when their model
# parameters and requested equivalent earths agree: they solve the same system.
function same_physical_state(z::EarthImpedance.Formula{:unified}, p::Formula{:unified})
    return same_physical_state(z.parameters, p.parameters) &&
           same_physical_state(z.equivalent_earth, p.equivalent_earth)
end

function initialize_buffers(
        selected::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        ::Type{T}, input, plan, buffers) where {T}
    buffers=initialize_buffers(selected.equivalent_earth, T, input, plan, buffers)
    buffers=initialize_buffers(SpectralIntegral, Val(:quad), T, input, plan, buffers)
    haskey(buffers, :current_map) && return buffers
    R=typeof(float(nominal(one(T))))
    n=length(plan.geometry.radius)
    axial_field=Matrix{Complex{T}}(undef, n, n)
    return merge(buffers,
        (
            axial_field, source_potential = similar(axial_field), current_map = similar(axial_field),
            enclosed_impedance = similar(axial_field), enclosed_potential = similar(axial_field),
            current_factor = similar(axial_field), current_rhs = similar(axial_field),
            radial_argument = Vector{Complex{T}}(undef, n), source_logscale = Vector{T}(undef, n),
            circumference_average = Vector{Complex{T}}(undef, n), radial_current = Vector{Complex{T}}(undef, n),
            earth_spectrum = (points = sizehint!(R[], 128), seeds = sizehint!(R[], 128),
                scales = sizehint!(R[], 16))))
end

# The Functor of the calculation at one frequency. Its state is the system that every pair
# shares, as plain values. The per-conductor arrays go to Unified's buffers, and the parts
# write the source coefficients into Unified's own two matrices.
function Functor(formula::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        input::NamedTuple; workspace)
    s=input.jω
    isfinite(s) && !iszero(s) || throw(DomainError(s, "jω must be finite and nonzero"))
    prescribed=formula.options.data.Γ
    longitudinal=prescribed isa Number ? prescribed : prescribed[input.frequency]
    for column in axes(input.rho, 2)
        validate(@view(input.rho[:, column]), formula,
            @view(input.epsilon[:, column]), @view(input.mu[:, column]),
            input.thickness)
    end
    for values in (input.rho, input.epsilon, input.mu)
        for column in axes(values, 2), row in axes(values, 1)

            same_physical_state(values[row, column], values[row, 1]) ||
                throw(ArgumentError("the unified earth-return calculation requires the same equivalent media for all conductor pairs"))
        end
    end
    Γ=oftype(s, longitudinal)
    sh=ntuple(m->conductivity(input.rho[m, 1])+s*input.epsilon[m, 1], 2)
    mu=ntuple(m->input.mu[m, 1], 2)
    gamma=ntuple(m->sqrt(s*mu[m]*sh[m]), 2)
    k2=ntuple(m->gamma[m]^2-Γ^2, 2)
    k=map(outgoing_root, k2)
    buffers=workspace.buffers
    geometry=workspace.plan.geometry
    for i in eachindex(buffers.radial_argument)
        medium=input.media[i]
        buffers.radial_argument[i]=k[medium]*geometry.radius[i]
        buffers.source_logscale[i]=abs(real(buffers.radial_argument[i]))
        buffers.circumference_average[i]=special_besselix(0, buffers.radial_argument[i])
        buffers.radial_current[i]=2 * (one(s)*π) * sh[medium] * geometry.radius[i]^2 *
                               bessel_current_ratio(buffers.radial_argument[i])
    end
    state=(jω = s, Γ, sh, mu, k2, k)
    destinations=(buffers.axial_field, buffers.source_potential)
    return Functor(formula, merge(input, (; destinations)), state)
end

"""
$(TYPEDSIGNATURES)

Convert the completed source coefficients to physical exterior impedance
\\[Ω/m\\] and potential coefficients \\[m/F\\] using the complete enclosed-current
map. With s=jω \\[1/s\\], prescribed Γ \\[1/m\\], and radial-current factors Fᵣ
\\[m/Ω\\],

```math
T_I=A_r^{-1}-F_r\\widetilde K,\\qquad
P_e T_I=\\widetilde P,\\qquad
Z_e T_I=\\widetilde K+\\Gamma^2\\widetilde P/s.
```

The work arrays contain a common exponential source-column scaling. Specifically, `current_map` stores TᵢD, not the unscaled Tᵢ. `circumference_average`
contains scaled I₀ values. The same D multiplies both source-coefficient arrays
and cancels from these right solves. The factorization is reused without
conjugating transposes. All conductors participate before selected physical
entries are copied by Engine.

This complete-field conversion includes both scalar- and vector-potential
contributions. The package then assembles cable contributions and
reductions before calculating total Y=jωP⁻¹ \\[S/m\\].

# Returns

- Named tuple of the physical `impedance` and `admittance` matrices: the impedance \\[Ω/m\\]
  and the potential coefficients \\[m/F\\].

# Reference

User-supplied manuscript, *Unified circumferentially averaged framework for
overhead, buried, and mixed conductor systems*, current and charge maps.
"""
function earth!(::Union{EarthImpedance.Formula{:unified}, Formula{:unified}}, functor::Functor,
        workspace)
    buffers=workspace.buffers
    u=functor.state
    for column in axes(buffers.axial_field, 2), row in axes(buffers.axial_field, 1)

        buffers.current_map[row, column]=(row==column ? inv(buffers.circumference_average[row]) :
                                       zero(u.jω))-
        buffers.radial_current[row]*buffers.axial_field[row, column]
    end
    copyto!(buffers.current_factor, transpose(buffers.current_map))
    factor=lu!(buffers.current_factor)
    copyto!(buffers.current_rhs, transpose(buffers.source_potential))
    ldiv!(factor, buffers.current_rhs)
    copyto!(buffers.enclosed_potential, transpose(buffers.current_rhs))
    @. buffers.enclosed_impedance=buffers.axial_field+u.Γ^2/u.jω*buffers.source_potential
    copyto!(buffers.current_rhs, transpose(buffers.enclosed_impedance))
    ldiv!(factor, buffers.current_rhs)
    copyto!(buffers.enclosed_impedance, transpose(buffers.current_rhs))
    return (impedance = buffers.enclosed_impedance, admittance = buffers.enclosed_potential)
end

:unified
