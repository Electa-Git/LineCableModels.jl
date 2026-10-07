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

function formulation_options(
        expression::Expression{<:Union{
            EarthImpedance.Formula{:unified}, Formula{:unified}}},
        ::Val{:Γ}, default, supplied)
    return validate(supplied, typeof(expression.selection))
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
    argument = validate(selected.options.data.Γ, typeof(selected))
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

Calculate the source-charge coefficient of conductor voltage \\[m/F\\] for
circumferentially averaged fields in two half-spaces. Air targets use the
interface z=0 as voltage reference. Buried targets use deep earth. With
source amplitudes q̃ and conductor voltages U,

```math
U=\\widetilde P\\widetilde q,\\qquad P_e T_I=\\widetilde P.
```

The coefficient includes the complete-field voltage relation, not just scalar
potential. Source-column exponential scaling matches the axial-field and
current-map arrays and cancels in the physical matrix solve. The air endpoint
and conductor field are combined before quadrature to preserve cancellation.

# Returns

- Scaled source-potential coefficient \\[m/F\\].

# Reference

User-supplied manuscript, *Unified circumferentially averaged framework for
overhead, buried, and mixed conductor systems*, voltage and source-charge maps.
"""
function source_potential_coefficient(::Union{Formula{:unified}, Val{:unified}},
        kind::Union{Val{:self}, Val{:mutual}}, source::Union{Val{1}, Val{2}},
        target::Union{Val{1}, Val{2}}, functor, pair, workspace)
    u=functor.state
    row, column=pair.row, pair.column
    hp, hq=abs(pair.heights[2]), abs(pair.heights[1])
    r=u.radius[row]
    average=u.circumference_average[row]
    sp, sq=u.source_logscale[row], u.source_logscale[column]
    integration=functor.options.data.integration
    context=(
        formula = :unified, frequency = imag(u.jω)/(2π), receiver = row, source = column)
    direct=earth_direct(
        kind, source, target, u, pair, r, average, u.radial_argument[row], sp, sq)
    return source_potential_coefficient(source, target, u, hp, hq, pair.separation,
        r, average, u.radial_argument[row], sp, sq,
        direct, integration, workspace.buffers; context)
end

function source_potential_coefficient(source::Val{S}, ::Val{2}, u, hp, hq, y,
        r, average, argument, sp, sq, direct, integration, numerical; context) where {S}
    value=u.jω/(one(u.jω)*π)*average*earth_spectral_term(Val(:voltage), Val(2), source,
        u, hp, hq, y, zero(r), sp+sq, integration.method,
        integration.options, numerical; context)
    return value+u.jω/(2*(one(u.jω)*π)*u.sh[2])*direct
end

function source_potential_coefficient(source::Val{S}, ::Val{1}, u, hp, hq, y,
        r, average, argument, sp, sq, direct, integration, numerical; context) where {S}
    πT=one(u.jω)*π
    if abs(real(nominal(argument)))<300
        R=typeof(float(nominal(real(u.jω))))
        padding=hq/2
        angle=min(R(π)/6, atan(R(nominal(hq))/(4max(R(nominal(y+r)), eps(R)))))
        angle=earth_contour_angle(u, angle)
        points=earth_spectral_points!(numerical.earth_spectrum, u, hq-padding, y, r, angle)
        g=(hp, hq, radius = r, padding, logscale = sq, i0minus = bessel_i0m1(argument))
        kernel=AirVoltageSpectrum{S, typeof(u), typeof(g)}(u, g)
        scale=max(R(abs(nominal(u.k[2]))), inv(R(nominal(hq))))
        contour=scale*cis(angle)
        integral=SpectralIntegral(t->contour*earth_spectral_value(
            kernel, contour*t, hq-padding, y, zero(r)))
        points ./= scale
        push!(points, one(scale))
        value,
        _=integrate(integral, integration.method, integration.options,
            numerical; points, coordinate_type = R,
            context = merge(context, (term = :air_voltage,)), observations = numerical.observations)
        return u.jω/πT*value+u.jω/(2πT*u.sh[1])*direct
    end
    value=u.jω/πT*average*earth_spectral_term(Val(:voltage), Val(1), source, u,
        hp, hq, y, zero(r), sp+sq, integration.method, integration.options, numerical; context)
    value+=u.jω/(2πT*u.sh[1])*direct
    return value+u.jω/πT*earth_spectral_term(Val(:air_reference), Val(1), source, u,
        zero(r), hq, y, r, sq, integration.method, integration.options, numerical; context)
end

function Expression(selected::Formula{:unified}, pair::EarthPair)
    return Expression(selected, source_potential_coefficient,
        Val(pair.row == pair.column ? :self : :mutual), Val.(layer_index(pair))...)
end

function formulation_options(::Expression{<:Formula{:unified},
        typeof(source_potential_coefficient),
        A}) where {
        A <: Tuple{
        Union{Val{:self}, Val{:mutual}}, Union{Val{1}, Val{2}}, Union{Val{1}, Val{2}}}}
    return FormulationOptions((Γ = 0, integration = (method = :quad, options = (;))))
end

function validate(reduction::EquivalentHomogeneous.Formula{:bottommost},
        ::Expression{<:Formula{:unified}, typeof(source_potential_coefficient)})
    return reduction
end

function earth_bindings(
        selected::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        model::EarthModel, reduction, physical::AbstractVector{<:EarthPair}, homogeneous,
        indices)
    binding=invoke(earth_bindings,
        Tuple{Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
            EarthModel, Any, AbstractVector{<:EarthPair}, Any, Any},
        selected, model, reduction, physical, homogeneous, collect(eachindex(physical)))
    options=first(binding.equations).declaration.options
    all(group->isequal(group.declaration.options, options), binding.equations) ||
        throw(ArgumentError("the unified earth-return calculation requires common formulation options for all conductor pairs"))
    groups=map(binding.equations) do group
        declaration=group.declaration
        primary=declaration.expression
        # Both coefficients are mathematical dependencies of either selected output.
        axial=primary.method === EarthImpedance.axial_field_coefficient ? primary :
              Expression(Val(:unified), EarthImpedance.axial_field_coefficient, primary.arguments...)
        potential=primary.method === source_potential_coefficient ? primary :
                  Expression(Val(:unified), source_potential_coefficient, primary.arguments...)
        (declaration = merge(declaration, (expression = (axial, potential),)),
            indices = group.indices)
    end
    return merge(binding, (equations = groups, output_indices = collect(indices)))
end

function earth_bindings(::Union{EarthImpedance.Formula{:unified}, Formula{:unified}},
        binding::NamedTuple, geometry::NamedTuple)
    inputs=map(binding.interactions) do interaction
        pair=interaction.pair
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
        (pair.row==pair.column, pair.layers, pair.heights, pair.separation,
            target_radius, source_radius)
    end
    return merge(binding, (reuse_inputs = inputs,))
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

# Formula-state construction consumes completed materials and engine-resolved layers.
function (selected::Union{EarthImpedance.Formula{:unified}, Formula{:unified}})(
        materials, binding, workspace, frequency::Int)
    s=workspace.input.jω[frequency]
    isfinite(s) && !iszero(s) || throw(DomainError(s, "jω must be finite and nonzero"))
    prescribed=first(binding.equations).declaration.options.data.Γ
    longitudinal=prescribed isa Number ? prescribed : prescribed[frequency]
    longitudinal isa Number && isfinite(longitudinal) ||
        throw(ArgumentError("Γ must be one finite scalar [1/m]"))
    for column in axes(materials.rho, 2)
        validate(@view(materials.rho[:, column]), selected,
            @view(materials.epsilon[:, column]), @view(materials.mu[:, column]),
            materials.thickness)
    end
    for values in (materials.rho, materials.epsilon, materials.mu)
        for column in axes(values, 2), row in axes(values, 1)

            same_physical_state(values[row, column], values[row, 1]) ||
                throw(ArgumentError("the unified earth-return calculation requires the same equivalent media for all conductor pairs"))
        end
    end
    Γ=oftype(s, longitudinal)
    sh=ntuple(m->conductivity(materials.rho[m, 1])+s*materials.epsilon[m, 1], 2)
    mu=ntuple(m->materials.mu[m, 1], 2)
    gamma=ntuple(m->sqrt(s*mu[m]*sh[m]), 2)
    k2=ntuple(m->gamma[m]^2-Γ^2, 2)
    k=map(outgoing_root, k2)
    buffers=workspace.buffers
    geometry=workspace.plan.geometry
    for i in eachindex(buffers.radial_argument)
        medium=binding.layers[i]
        buffers.radial_argument[i]=k[medium]*geometry.radius[i]
        buffers.source_logscale[i]=abs(real(buffers.radial_argument[i]))
        buffers.circumference_average[i]=special_besselix(0, buffers.radial_argument[i])
        buffers.radial_current[i]=2 * (one(s)*π) * sh[medium] * geometry.radius[i]^2 *
                               bessel_current_ratio(buffers.radial_argument[i])
    end
    state=(jω = s, Γ, sh, mu, k2, k, radius = geometry.radius,
        radial_argument = buffers.radial_argument, source_logscale = buffers.source_logscale,
        circumference_average = buffers.circumference_average)
    return (coefficients = (buffers.axial_field, buffers.source_potential), state)
end

function (selected::EarthImpedance.Formula{:unified})(state::NamedTuple, interaction::NamedTuple, declaration)
    binding=(pair = interaction.pair, physical_pair = interaction.physical_pair,
        kind = declaration.kind, expression = declaration.expression)
    return EarthImpedance.Functor(binding, state, declaration.options)
end
function (selected::Formula{:unified})(state::NamedTuple, interaction::NamedTuple, declaration)
    binding=(pair = interaction.pair, physical_pair = interaction.physical_pair,
        kind = declaration.kind, expression = declaration.expression)
    return Functor(binding, state, declaration.options)
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

- Named tuple of physical `impedance` and `potential` matrices.

# Reference

User-supplied manuscript, *Unified circumferentially averaged framework for
overhead, buried, and mixed conductor systems*, current and charge maps.
"""
function earth!(::Union{EarthImpedance.Formula{:unified}, Formula{:unified}}, calculation, workspace)
    buffers=workspace.buffers
    u=calculation.state
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
    return (impedance = buffers.enclosed_impedance, potential = buffers.enclosed_potential)
end

function earth_bindings(z::EarthImpedance.Formula{:unified}, p::Formula{:unified}, impedance, admittance)
    same_physical_state(z.parameters, p.parameters) &&
    same_physical_state(z.equivalent_earth, p.equivalent_earth) &&
    same_physical_state(first(impedance.equations).declaration.options.data,
        first(admittance.equations).declaration.options.data) || return nothing
    return (; impedance, admittance)
end

:unified
