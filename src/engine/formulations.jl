# Engine-owned formulation hierarchy.
"""
$(TYPEDEF)

Select the LineCableModels backend for concentric coaxial cable assemblies.

Series impedance uses the equivalent concentric representation supplied by
DataModel. The default local shunt model uses the equivalent annular layer. Explicit
`shunt_model=:boundary` resolves eligible open wires and finite tapes in a
lossless, radially layered circular shielded domain during blueprint
construction. Other geometry and material selections retain their
equivalent-coaxial treatment.
"""
struct LineCableModelsCoaxial end

"""Identify the coaxial backend without executing or configuring a computation."""
description(::Type{LineCableModelsCoaxial}; compact::Bool = false) = "coaxial"
function description(::LineCableModelsCoaxial; compact::Bool = false)
    description(LineCableModelsCoaxial; compact)
end
formula_id(::Type{LineCableModelsCoaxial}) = :coaxial
formula_id(::LineCableModelsCoaxial) = :coaxial

"""
$(TYPEDEF)

Select the Julia-native Gmsh/GetDP finite-element backend.

`options` stores the field-model choice and matrix reductions. Set
`options=(physics=:quasi_tem,)` (default) or `(physics=:quasi_fw,)`.
Execution controls belong to `compute(...; options=(...))` and are validated
by [`computation_options`](@ref) for `LineCableModelsFEM`.

The quasi-TEM model solves independent axial ``A_z/u`` and scalar electric
Helmholtz blocks in one factorization. The magnetic excitation is one ampere.
The electric excitation is one ampere per meter. Both models retain conduction
and displacement through ``κ=σ+jωε`` \\[S/m\\], where ``ω`` is angular
frequency \\[rad/s\\], ``σ`` conductivity \\[S/m\\], and ``ε`` permittivity \\[F/m\\].

The quasi-full-wave model uses phasors ``e^{jωt-Γz}`` and expands
``A_z=a``, ``A_t=Γb``, ``φ=Γv`` before taking ``Γ→0``. Here ``Γ`` is
the longitudinal propagation constant \\[1/m\\], ``a`` has units \\[T m\\],
``b`` \\[T m²\\], and ``v`` \\[V m\\]. In the media outside the equipotential
metal terminals, the retained equations are

```math
-\\nabla_t\\!\\cdot(\\nu\\nabla_t a)+j\\omega\\kappa a=0,
\\qquad
C^*(\\nu Cb)+j\\omega\\kappa b+\\kappa\\nabla_t v-\\nu\\nabla_t a=0,
\\qquad
-\\nabla_t\\!\\cdot[\\kappa(j\\omega b+\\nabla_t v)]+j\\omega\\kappa a=0.
```

Here ``\\nu=1/\\mu`` \\[m/H\\], ``μ`` is permeability \\[H/m\\],
``Cb=∂_x b_y-∂_y b_x`` and ``C^*h=(∂_y h,-∂_x h)``.
The finite-conductivity axial ``a/u`` block supplies both the series voltage
drop and the normalized transverse-current source. There is no independent
electric excitation. A tree gauge removes the gradient freedom of ``b``.
Its tangential trace vanishes on every terminal contour and the entire outer
boundary. ``v=0`` on the earth-side outer reference, with natural electric
conditions on the air side.

Terms of order ``Γ^2`` in the axial equation are omitted. The transverse
Ampère equation remains at the order used to extract shunt response: continuity
alone constrains a divergence and cannot supply that vector balance across a
material interface. Choosing a smaller numerical ``Γ`` cannot restore it.

Terminal voltage includes the vector potential. For a path ``ℓ_i`` oriented
from the earth reference to terminal ``i``, the inverse-admittance matrix is

```math
P_{ij}=\\frac{v_i-v_{\\rm ref}+j\\omega\\int_{\\ell_i}b\\cdot d\\ell}{I_j},
\\qquad Y=P^{-1},\\qquad P_e=j\\omega P.
```

``I_j`` is the imposed axial current \\[A\\], ``P`` has units \\[Ω m\\],
``Y`` \\[S/m\\], and the analytical potential coefficient ``P_e`` \\[m/F\\].
The backend uses physical vertical paths to the lowest mesh node of each
terminal contour, with zero transverse field inside equipotential metal.
It integrates the pulled-back edge field through the infinite shell.
The scalar trace ``v_i/I_j`` alone is gauge dependent and is saved separately
as `Pscalar.tsv` diagnostics. Series impedance is ``Z_{ij}=-U_i/I_j`` \\[Ω/m\\],
where ``U_i`` is the axial electric unknown \\[V/m\\].

This is a two-dimensional first-order reduction of Maxwell's potential
equations, retaining transverse induction and displacement. It does not solve
for a finite propagation constant or constitute the Darwin approximation.
The full-vector potential equations and the role of terminal conditions and
gauging are described by G. Ciuprina and R. V. Sabriego, *Electric circuit
element boundary conditions for electromagneto-quasistatic and full wave
models in A, φ potentials and their finite element implementation*, Journal
of Mathematics in Industry **14**, 27 (2024),
[doi:10.1186/s13362-024-00165-6](https://doi.org/10.1186/s13362-024-00165-6),
Sect. 4 and Appendix B. The longitudinal reduction and vertical voltage-path
convention above are specific to this backend. The paper validates 3D models.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEM{M <: NamedTuple, O <: FormulationOptions, D <: NamedTuple} <:
       AbstractFormulation
    "Shared scientific formula selections, independent of FEM execution controls."
    methods::M
    "Field model and line-parameter matrix reductions."
    options::O
    "Requested formula definitions."
    definitions::D
end

"""Identify the FEM backend without loading meshing or solver packages."""
description(::Type{<:LineCableModelsFEM}; compact::Bool = false) = "FEM"
function description(::LineCableModelsFEM; compact::Bool = false)
    description(LineCableModelsFEM; compact)
end
formula_id(::Type{<:LineCableModelsFEM}) = :fem
formula_id(::LineCableModelsFEM) = :fem
formulation_options(value::LineCableModelsFEM) = value.options
function Base.pairs(value::LineCableModelsFEM; quantity = nothing)
    pairs(LineCableModelsFEM,
        (methods = value.methods,
            requested = value.definitions,
            options = value.options.data);
        quantity)
end

"""FEM's coupled field equations retain all four constitutive selections."""
function Base.pairs(::Type{LineCableModelsFEM}; quantity = nothing)
    return pairs((insulation_admittance = InsulationAdmittance.Formula,
        semicon_admittance = SemiconAdmittance.Formula,
        earth_properties = Earth.FrequencyDependent.Formula,
        temperature_dependence = TemperatureDependent.Formula))
end
function description(::Type{LineCableModelsFEM}, slot::Val; compact::Bool=false)
    description(LineParametersFormulation, slot; compact)
end
function Base.pairs(::Type{LineCableModelsFEM}, retained::NamedTuple; quantity = nothing)
    pairs(LineParametersFormulation, retained; quantity, owner = LineCableModelsFEM)
end

"""
$(TYPEDEF)

Report a failure in finite-element model construction, meshing, solving, or result validation.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEMError <: Exception
    "Failure category."
    category::Symbol
    "Stable identifier of the object associated with the failure."
    object_id::String
    "Field or derived datum that failed validation."
    field::Symbol
    "Human-readable failure description."
    message::String
    "Retained run directory, or `nothing` before run creation."
    run_directory::Union{Nothing, String}
end

function LineCableModelsFEMError(
        category::Symbol,
        object_id,
        field::Symbol,
        message::AbstractString;
        run_directory::Union{Nothing, AbstractString} = nothing
)
    path = run_directory === nothing ? nothing : String(run_directory)
    return LineCableModelsFEMError(
        category, String(object_id), field, String(message), path
    )
end

function Base.showerror(io::IO, error::LineCableModelsFEMError)
    print(
        io,
        "LineCableModelsFEMError(",
        error.category,
        ", object=",
        repr(error.object_id),
        ", field=:",
        error.field,
        "): ",
        error.message
    )
    error.run_directory === nothing || print(
        io, "; retained run directory: ", error.run_directory
    )
end

"""
$(TYPEDEF)

Supertype for Engine impedance formulations.
"""
abstract type AbstractImpedanceFormulation <: AbstractFormulation end
"""
Select equations for the surface impedance of conductors [Ω/m]. Concrete subtypes build, in a
`Functor` method, the values that inner, outer and transfer surface impedances share, from
the conductor dimensions and material properties. They implement
`InternalImpedance.internal_impedance` for their supported `inner`, `outer`, and `transfer`
cases.
"""
abstract type InternalImpedanceFormulation <: AbstractImpedanceFormulation end
"""
Select the pipe contribution and its backend and topology applicability. Concrete
subtypes extend `validate(design, selected, backend)`, which admits the topology of
`design` or throws. An impedance equation must be supplied separately from admission.
"""
abstract type PipeImpedanceFormulation <: AbstractImpedanceFormulation end
"""
Select magnetic impedance across an insulation annulus [Ω/m]. Concrete subtypes
implement `InsulationImpedance.insulation_impedance` on their selected type.
"""
abstract type InsulationImpedanceFormulation <: AbstractImpedanceFormulation end
"""
Select earth-return impedance equations [Ω/m]. Concrete subtypes implement
indexed `EarthImpedance.earth_impedance` methods on evaluated material inputs.
Air, earth, and mixed selections refer to conductor locations.
"""
abstract type EarthImpedanceFormulation <: AbstractImpedanceFormulation end

"""Supertype for local shunt geometry and dielectric and earth admittance selections."""
abstract type AbstractAdmittanceFormulation <: AbstractFormulation end
"""
Select cable-local shunt geometry, independently of material admittivity.
Concrete subtypes implement `internal_shunt_response` during blueprint construction.
"""
abstract type ShuntModelFormulation <: AbstractAdmittanceFormulation end
"""
Select insulation admittivity [S/m]. Concrete subtypes implement
`InsulationAdmittance.insulation_material`. Radial geometry is applied separately.
"""
abstract type InsulationAdmittanceFormulation <: AbstractAdmittanceFormulation end
"""
Select semiconducting-layer admittivity [S/m]. Concrete subtypes implement
`SemiconAdmittance.semicon_material`. Radial geometry is applied separately.
"""
abstract type SemiconAdmittanceFormulation <: AbstractAdmittanceFormulation end
"""
Select earth potential-coefficient equations [m/F]. Concrete subtypes implement
indexed `EarthAdmittance.earth_potential_coefficient` methods on evaluated
material inputs. Matrix assembly converts
the potential coefficients to shunt admittance [S/m].
"""
abstract type EarthAdmittanceFormulation <: AbstractAdmittanceFormulation end

"""
Validate the admittivity [S/m] that a dielectric law returns, before radial aggregation
represents it as `Complex{T}`. A result requiring a wider scalar type than `T` is rejected
instead of silently discarding precision or uncertainty. Return `value`.
"""
function validate(value, ::Union{InsulationAdmittanceFormulation, SemiconAdmittanceFormulation},
        ::Type{T}) where {T <: Real}
    value isa Number && !(value isa Bool) && isfinite(value) || throw(DomainError(
        value, "a dielectric material law must return a finite scalar admittivity [S/m]"))
    promote_type(T, typeof(real(value))) === T || throw(ArgumentError(
        "dielectric admittivity requires scalar type $(typeof(real(value))); " *
        "use a material and problem scalar type that preserves its precision and uncertainty"))
    return value
end

# The expression of an earth formula for one interaction: its kind and its source and target
# layers.
function Expression(formula::EarthImpedanceFormulation, pair::EarthPair)
    return Expression(formula, EarthImpedance.earth_impedance,
        Val(pair.row == pair.column ? :self : :mutual), Val.(layer_index(pair))...)
end
function Expression(formula::EarthAdmittanceFormulation, pair::EarthPair)
    return Expression(formula, EarthAdmittance.earth_potential_coefficient,
        Val(pair.row == pair.column ? :self : :mutual), Val.(layer_index(pair))...)
end

"""
$(TYPEDSIGNATURES)

Build the expression of an earth recipe for `pair`: the `air` formula for layers (1, 1), the
`earth` formula for (2, 2) and the `mixed` formula for (1, 2) and (2, 1). The pair's layers
are the decided ones: a recipe describes a homogeneous earth.

# Errors

- Throws `ArgumentError` when the recipe has no formula for the pair's route, or for another
  layer pair.
"""
function Expression(recipe::NamedTuple, pair::EarthPair)
    source, target = layer_index(pair)
    route = (source, target) == (1, 1) ? :air : (source, target) == (2, 2) ? :earth :
            (source, target) in ((1, 2), (2, 1)) ? :mixed :
            throw(ArgumentError("homogeneous selection is not defined for source in layer " *
                "$source and target in layer $target"))
    selected = get(recipe, route, nothing)
    selected === nothing &&
        throw(ArgumentError("explicit equation recipe has no requested :$route case"))
    return Expression(selected, pair)
end

# Equation-specific geometric restrictions extend the existing validation protocol.
validate(pair::EarthPair, ::Expression) = pair

# Whether the formula of `expression` admits earth layer `k`. It does when a method of its
# operation accepts `Val{k}` in the source or the target position, with the `Functor` and the
# workspace that the evaluation passes, whatever the types of the other arguments. A layer
# left generic, such as `::Val{S}`, admits every layer. The earth decision, the error message
# of `validate(expression, layers)` and the standalone formula call read the layers here.
function _admits_layer(expression::Expression, k::Int)
    F = typeof(expression.selection)
    return !isempty(methods(expression.method, Tuple{F, Any, Val{k}, Any, Functor, Any})) ||
           !isempty(methods(expression.method, Tuple{F, Any, Any, Val{k}, Functor, Any}))
end

"""
$(TYPEDSIGNATURES)

Check that the formula of `expression` defines it for an earth with one entry of `layers` per
medium, before evaluating it: the layers of an earth model, or the resistivities of the
standalone formula call. The signatures of a formula's methods declare the earth layers it
handles, whatever the types of their runtime arguments. In layers that the formula admits,
`validate(expression)` checks the expression. Return `expression`.

# Errors

- Throws `ArgumentError` with the layer count N and the highest layer up to N that the
  formula admits, when the formula does not admit the source or target layer.
- Throws the `ArgumentError` of `validate(expression)` for a missing expression in layers
  that the formula admits.

When the formula admits layer N + 1, the message reads "defined beyond layer N" instead.
"""
function validate(expression::Expression{<:Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}},
        layers::Union{Tuple, AbstractVector})
    selected = expression.selection
    signature = Tuple{typeof(selected), map(typeof, expression.arguments)..., Functor, Any}
    hasmethod(expression.method, signature) && return expression
    interaction(::Val{K}, ::Val{S}, ::Val{T}) where {K, S, T} = (K, S, T)
    kind, source, target = interaction(expression.arguments...)
    if _admits_layer(expression, source) && _admits_layer(expression, target)
        validate(expression)
        return expression
    end
    N = length(layers)
    defined = _admits_layer(expression, N + 1) ? "defined beyond layer $N" :
              "defined up to layer $(something(findlast(k -> _admits_layer(expression, k), 1:N), 0))"
    throw(ArgumentError(
        "the earth model has $N layers and formula " *
        ":$(formula_id(selected)) is $defined; it has no expression " *
        "for a $kind interaction from layer $source to layer $target"))
end

function validate(earth::EarthModel, formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation})
    validate(earth)
    earth.vertical_layers &&
        throw(ArgumentError("earth-return equations require horizontal interfaces or an explicit EquivalentHomogeneous reduction"))
    return earth
end

function validate(rho::AbstractVector,
        formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        epsilon::AbstractVector, mu::AbstractVector, thickness)
    length(rho) == length(epsilon) == length(mu) ||
        throw(DimensionMismatch("material vectors must align"))
    all(x -> x > 0 && !isnan(x), rho) ||
        throw(DomainError(rho, "resistivities must be positive, including infinite air resistivity"))
    all(x -> isfinite(x) && !iszero(x), epsilon) && all(x -> isfinite(x) && x > 0, mu) ||
        throw(DomainError((epsilon, mu),
            "permittivities must be nonzero and finite; permeabilities positive and finite"))
    epsilon[1] > 0 || throw(DomainError(epsilon[1], "air permittivity must be positive"))
    all(>(0), epsilon) ||
        throw(DomainError(epsilon,
            "formula :$(formula_id(formula)) requires positive permittivity; the earth-material data type permits artificial negative values"))
    if thickness === nothing
        length(rho) == 2 || throw(DimensionMismatch(
            "material vectors without layer thicknesses describe air and one earth medium"))
    else
        length(thickness) == length(rho) ||
            throw(DimensionMismatch("layer thicknesses must align with materials"))
        isinf(first(thickness)) && isinf(last(thickness)) &&
        all(x -> isfinite(x) && x > 0, @view(thickness[2:(end - 1)])) ||
            throw(DomainError(thickness,
                "air and bottom half-spaces must be infinite; internal layers positive and finite"))
    end
    return rho
end

"""
$(TYPEDSIGNATURES)

Return the normalized options of the surface impedance `kind` of the internal-impedance
formula `formula`: the section of its options for this kind, or empty options.
"""
function formulation_options(formula::InternalImpedanceFormulation, ::Val{Kind}) where {Kind}
    return FormulationOptions(get(formula.options.data, Kind, (;)))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of an insulation-impedance formula for one insulation annulus at one
frequency, after checking its radii, its relative permeability and `jω`. The state is empty.
"""
function Functor(formula::InsulationImpedanceFormulation, input::NamedTuple;
        workspace = nothing)
    (; r_in, r_ex, mu_r, jω) = input
    T = typeof(r_in)
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in <= r_ex ||
        throw(DomainError((r_in, r_ex), "insulation radii must satisfy 0 ≤ r_in ≤ r_ex [m]"))
    isfinite(mu_r) && mu_r > zero(T) || throw(DomainError(mu_r,
        "relative insulation permeability must be positive and finite"))
    isfinite(jω) || throw(DomainError(jω, "jω must be finite"))
    return Functor(formula, input, (;))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of an insulation admittivity law for one material at one frequency and
temperature, after checking both. The state is empty.
"""
function Functor(formula::InsulationAdmittanceFormulation, input::NamedTuple;
        workspace = nothing)
    (; frequency, temperature) = input
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency, "insulation constitutive frequency must be positive and finite"))
    isfinite(temperature) || throw(DomainError(
        temperature, "insulation constitutive temperature must be finite"))
    return Functor(formula, input, (;))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of a semiconductor admittivity law for one material at one frequency and
temperature, after checking both. The state is empty.
"""
function Functor(formula::SemiconAdmittanceFormulation, input::NamedTuple;
        workspace = nothing)
    (; frequency, temperature) = input
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency, "semicon constitutive frequency must be positive and finite"))
    isfinite(temperature) || throw(DomainError(
        temperature, "semicon constitutive temperature must be finite"))
    return Functor(formula, input, (;))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of an earth calculation of the impedance formula `formula` at one
frequency, whose parts write the impedance coefficients of each conductor pair into `Zearth`
of the workspace buffers, with an empty state.
"""
function Functor(formula::EarthImpedanceFormulation, input::NamedTuple; workspace)
    return Functor(formula, merge(input, (destinations = (workspace.buffers.Zearth,),)), (;))
end

"""
$(TYPEDSIGNATURES)

Build the Functor of an earth calculation of the admittance formula `formula` at one
frequency, whose parts write the potential coefficients of each conductor pair into `Pearth`
of the workspace buffers, with an empty state.
"""
function Functor(formula::EarthAdmittanceFormulation, input::NamedTuple; workspace)
    return Functor(formula, merge(input, (destinations = (workspace.buffers.Pearth,),)), (;))
end

"""
$(TYPEDSIGNATURES)

Check the input of one conductor pair before an earth formula evaluates it: a finite nonzero
`jω`, the aligned material vectors of the pair, and, on a layered earth, the pair's heights
against the layer thicknesses. Return `input`.
"""
function validate(input::NamedTuple,
        formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation})
    isfinite(input.jω) && !iszero(input.jω) ||
        throw(DomainError(input.jω, "jω must be finite and nonzero"))
    validate(input.rho, formula, input.epsilon, input.mu, input.thickness)
    input.thickness === nothing || validate(input.pair, input.thickness)
    return input
end

"""
$(TYPEDSIGNATURES)

Evaluate the earth coefficient of one indexed conductor pair at fixed angular frequency.
Material vectors list air and every physical earth layer when layer thicknesses are given,
and exactly `(air,soil)` otherwise. The call runs the checks of a computation: the pair and
its geometry, the expression for its kind and layers on an earth with these media, the
formula's options for that expression, and the pair's input. Return the coefficient as
`Complex{T}`.

# Arguments

- `resistivity`: aligned resistivities [Ω·m].
- `permittivity`: absolute permittivities [F/m].
- `permeability`: absolute permeabilities [H/m].
- `jω`: imaginary angular frequency [1/s].
- `pair`: indexed conductor interaction. Lengths [m].

# Keywords

- `thickness`: aligned layer thicknesses [m] for a layered earth.
- `physical_pair`: the physical pair before an equivalent-earth reduction.
- `workspace`: the computation workspace, when the expression reads its buffers.
"""
function (formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation})(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair;
        thickness = nothing, physical_pair = pair, workspace = nothing
) where {T <: Real}
    validate(pair)
    expression = Expression(formula, pair)
    validate(pair, expression)
    validate(expression, resistivity)
    options = only(formulation_options(formula, (expression,)).options)
    # One pair does not need destinations. The plain constructor gives the empty state of a
    # formula that does not share values. Unified, whose state is its whole system, runs only in
    # a computation.
    functor = Functor(formula, (; jω, thickness, pair, physical = physical_pair,
        rho = resistivity, epsilon = permittivity, mu = permeability, options), (;))
    validate(functor.input, formula)
    value = expression(functor, workspace)
    value isa Number && isfinite(value) ||
        throw(DomainError((value,), "earth coefficients must be finite scalars"))
    return oftype(jω, value)
end

"""
$(TYPEDSIGNATURES)

Build the backend formulation whose identity is `tag`: `:coaxial`, `:cable_constants`,
`:fem` or `:pscad`. `Formulation(::Val{tag}; kwargs...)` forwards `kwargs` to the
backend's constructor. PSCAD defines its own method.

# Errors

- Throws `ArgumentError` for an unknown `tag`.
"""
Formulation(tag::Symbol; kwargs...) = Formulation(Val(tag); kwargs...)

Formulation(::Val{tag}; kwargs...) where {tag} =
    throw(ArgumentError("unknown backend formulation :$tag"))

function _fem_formulation(
        insulation_admittance, semicon_admittance, earth_properties, temperature_dependence,
        options::FormulationOptions
)
    methods = (
        insulation_admittance = InsulationAdmittance.Formula(insulation_admittance),
        semicon_admittance = SemiconAdmittance.Formula(semicon_admittance),
        earth_properties = earth_properties === nothing ? nothing :
                           Earth.FrequencyDependent.Formula(earth_properties),
        temperature_dependence = temperature_dependence === nothing ? nothing :
                                 TemperatureDependent.Formula(temperature_dependence)
    )
    definitions = (; insulation_admittance, semicon_admittance, earth_properties,
        temperature_dependence)
    return LineCableModelsFEM(methods, formulation_options(LineCableModelsFEM, options),
        definitions)
end

"""
$(TYPEDSIGNATURES)

Construct the Gmsh/GetDP finite-element formulation. FEM defines its field equations
and selects four material laws. Each law and `options` accepts a
scalar or an explicit `Grid`/`Gridspace`. Varying inputs return a
`Gridspace{LineCableModelsFEM}`.

# Keywords

- `insulation_admittance`: insulation admittivity law. `:default` routes to
  `:lossless`.
- `semicon_admittance`: semicon admittivity law. `:default` routes to
  `:lossless`.
- `earth_properties`: soil frequency-dependent constitutive law. `:default`
  routes to the explicit `:constant` pass-through, while `nothing` preserves
  the declared static soil. Equivalent-earth reductions are unsupported. Air
  uses its declared static properties.
- `temperature_dependence`: cable-material resistivity law. `:default` selects
  the linear law and `nothing` retains reference resistivity. Operating
  temperature belongs to `LineParametersProblem`.
- `options=(;)`: field model (`physics=:quasi_tem` or `:quasi_fw`) and bundle,
  Kron, and ideal-transposition reductions. Hyphenated strings and symbols
  are also accepted for `physics`. Julia parses `:quasi-fw` as subtraction.
  use `:quasi_fw` or `Symbol("quasi-fw")`.
  Pass execution controls to `compute(...; options=(...))`.
- `combine=:product`: product or zip composition among varying inputs.

Analytical impedance and admittance kernel keywords are rejected. Supported enclosure
geometry is represented directly in the FEM domain.
"""
function LineCableModelsFEM(;
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
        earth_properties = formula(:default),
        temperature_dependence = formula(:default),
        options = FormulationOptions(),
        combine::Symbol = :product
)
    return parameterize(
        LineCableModelsFEM,
        (inputs...) -> _fem_formulation(inputs[1:(end - 1)]...,
            last(inputs) isa NamedTuple ? FormulationOptions(last(inputs)) : last(inputs)),
        (insulation_admittance, semicon_admittance, earth_properties,
            temperature_dependence, options);
        combine
    )
end

Formulation(::Val{:fem}; kwargs...) = LineCableModelsFEM(; kwargs...)

# An earth equation admits an equivalent-earth reduction only through its own method.
function validate(reduction::EquivalentHomogeneous.AbstractRule,
        expression::Expression{<:Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}})
    throw(ArgumentError("$expression does not admit equivalent-earth reduction :$(formula_id(reduction))"))
end

"""Expose FEM constitutive and admittance selections, field model, and reductions."""
function Base.NamedTuple(value::LineCableModelsFEM)
    record = function (selected)
        selected === nothing && return nothing
        selected isa Symbol && return NamedTuple(formula(selected))
        selected isa NamedTuple && return map(record, selected)
        return NamedTuple(selected)
    end
    Record=NamedTuple{(:backend, :requested, :methods, :options),
        Tuple{Symbol, NamedTuple, NamedTuple, NamedTuple}}
    return Record((
        :fem, map(record, value.definitions), map(record, value.methods), value.options.data))
end
