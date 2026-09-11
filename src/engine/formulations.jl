# Engine-owned formulation hierarchy.
"""
$(TYPEDEF)

Select the LineCableModels backend for concentric coaxial cable assemblies.

Nonconcentric cable parts must reach this backend through an equivalent
concentric representation supplied by DataModel.
"""
struct LineCableModelsCoaxial end

"""
$(TYPEDEF)

Select the Julia-native Gmsh/GetDP finite-element backend.

`options` stores the field-model choice and matrix reductions. Set
`options=(physics=:quasi_tem,)` (default) or `(physics=:quasi_fw,)`.
Execution controls belong to `compute(...; options=(...))` and are validated
by [`computation_options`](@ref) for `LineCableModelsFEM`.

The quasi-TEM model solves independent axial ``A_z/u`` and scalar electric
Helmholtz blocks in one factorization. The magnetic excitation is one ampere;
the electric excitation is one ampere per metre. Both models retain conduction
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
drop and the normalized transverse-current source; there is no independent
electric excitation. A tree gauge removes the gradient freedom of ``b``.
Its tangential trace vanishes on every terminal contour and the entire outer
boundary; ``v=0`` on the earth-side outer reference, with natural electric
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
convention above are specific to this backend; the paper validates 3D models.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEM{M <: NamedTuple, O <: FormulationOptions, D <: NamedTuple} <:
       AbstractFormulation
    "Shared scientific formula selections, independent of FEM execution controls."
    methods::M
    "Field model and line-parameter matrix reductions."
    options::O
    "Requested formula definitions retained for provenance."
    definitions::D
end

"""
$(TYPEDEF)

Report a finite-element adaptation, mesh, solve, or result-contract failure.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEMError <: Exception
    "Failure category."
    category::Symbol
    "Stable identifier of the object that owns the failure."
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
abstract type InternalImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type PipeImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type InsulationImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type EarthImpedanceFormulation <: AbstractImpedanceFormulation end

abstract type AbstractAdmittanceFormulation <: AbstractFormulation end
abstract type InsulationAdmittanceFormulation <: AbstractAdmittanceFormulation end
abstract type SemiconAdmittanceFormulation <: AbstractAdmittanceFormulation end
abstract type EarthAdmittanceFormulation <: AbstractAdmittanceFormulation end

"Return whether an earth formulation consumes homogeneous or stratified media."
function media end

"Declare source-owned physical hook defaults and admitted overrides for an equation binding."
function hooks end

"Resolve scalar or homogeneous three-field selections through their formula owner."
function Formulation(::Type{F}, selected) where {F <: Union{
        EarthImpedanceFormulation, EarthAdmittanceFormulation}}
    return F(selected)
end

function Formulation(::Type{F}, selected::NamedTuple) where {F <: Union{
        EarthImpedanceFormulation, EarthAdmittanceFormulation}}
    names = (:air, :earth, :mixed)
    length(selected) == 3 && all(in(keys(selected)), names) || throw(ArgumentError(
        "homogeneous earth selections require exactly air, earth and mixed"))
    return map(F, NamedTuple{names}(selected))
end

"Resolve the selected formula for exact source and target layer indices."
Formulation(selected::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
    ::Val{S}, ::Val{T}) where {S, T} = selected

Formulation(selected::NamedTuple{(:air, :earth, :mixed)}, ::Val{1}, ::Val{1}) = selected.air
Formulation(selected::NamedTuple{(:air, :earth, :mixed)}, ::Val{2}, ::Val{2}) = selected.earth
Formulation(selected::NamedTuple{(:air, :earth, :mixed)}, ::Val{1}, ::Val{2}) = selected.mixed
Formulation(selected::NamedTuple{(:air, :earth, :mixed)}, ::Val{2}, ::Val{1}) = selected.mixed

function Formulation(::NamedTuple{(:air, :earth, :mixed)}, ::Val{S}, ::Val{T}) where {S, T}
    throw(ArgumentError(
        "homogeneous selection is not defined for source in layer $S and target in layer $T"))
end

function validate(selected::NamedTuple{(:air, :earth, :mixed)}, earth::EarthModel)
    validate(earth)
    !earth.vertical_layers && length(earth.layers) == 2 &&
        all(layer -> isinf(layer.thickness), earth.layers) || throw(ArgumentError(
        "air/earth/mixed selections require physical air and one homogeneous soil half-space; use a scalar formulation for a layered model"))
    return selected
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        pair::EarthPair)
    return only(validate(formula, (pair,)))
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        pairs::Union{Tuple, AbstractVector{<:EarthPair}})
    equations = map(pairs) do pair
        validate(pair)
        equation = validate(FormulaMethod(formula, pair))
        validate(pair, equation)
        equation
    end
    identities = unique(equations)
    bindings = map(identities) do equation
        declared = hooks(equation)
        all(in(declared.configurable), keys(formula.hooks)) || throw(ArgumentError(
            "an explicit physical hook is unused by $equation"))
        selected_hooks = merge(declared.defaults, formula.hooks)
        defaults = selected_hooks.contribution === nothing ? computation_options(equation) :
                   computation_options(equation, selected_hooks.contribution)
        (equation = equation, kind = typeof(first(equation.arguments)).parameters[1],
            hooks = selected_hooks, defaults = defaults)
    end
    admitted = union((keys(binding.defaults) for binding in bindings)...)
    unknown = setdiff(keys(formula.options), admitted)
    isempty(unknown) || throw(ArgumentError(
        "unused numerical sections $(Tuple(unknown)) for required cases of :$(formula_id(formula))"))
    resolved = map(bindings) do binding
        names = Tuple(intersect(keys(formula.options), keys(binding.defaults)))
        options = computation_options(binding.equation, binding.defaults, formula.options[names])
        (equation = binding.equation, kind = binding.kind,
            hooks = binding.hooks, options = options)
    end
    return map(equation -> resolved[findfirst(==(equation), identities)], equations)
end

# Equation-specific geometric restrictions extend the existing validation protocol.
validate(pair::EarthPair, ::FormulaMethod) = pair

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}, count::Integer)
    count in formula.assumptions.layers || throw(DimensionMismatch(
        "formula :$(formula_id(formula)) requires $(formula.assumptions.layers) media including air; received $count"))
    return formula
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}, earth::EarthModel)
    validate(earth)
    earth.vertical_layers &&
        throw(ArgumentError("earth-return equations require horizontal interfaces or an explicit EquivalentHomogeneous reduction"))
    validate(formula, length(earth.layers))
    return formula
end

function validate(formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        rho::AbstractVector, epsilon::AbstractVector, mu::AbstractVector, thickness)
    length(rho) == length(epsilon) == length(mu) ||
        throw(DimensionMismatch("material vectors must align"))
    validate(formula, length(rho))
    all(x -> x > 0 && !isnan(x), rho) ||
        throw(DomainError(rho, "resistivities must be positive, including infinite air resistivity"))
    all(x -> isfinite(x) && !iszero(x), epsilon) && all(x -> isfinite(x) && x > 0, mu) ||
        throw(DomainError((epsilon, mu),
            "permittivities must be nonzero and finite; permeabilities positive and finite"))
    epsilon[1] > 0 || throw(DomainError(epsilon[1], "air permittivity must be positive"))
    restriction = formula.assumptions.permittivity
    restriction in (:positive, :nonzero) ||
        throw(ArgumentError("unknown source permittivity restriction"))
    restriction === :positive && !all(>(0), epsilon) &&
        throw(DomainError(epsilon,
            "formula :$(formula_id(formula)) requires positive permittivity; the earth-material data type permits artificial negative values"))
    if thickness === nothing
        media(formula) === Val(:stratified) && length(rho) > 2 &&
            throw(DimensionMismatch("stratified equations require aligned physical layer thicknesses"))
    else
        length(thickness) == length(rho) ||
            throw(DimensionMismatch("layer thicknesses must align with materials"))
        isinf(first(thickness)) && isinf(last(thickness)) &&
        all(x -> isfinite(x) && x > 0, @view(thickness[2:(end - 1)])) ||
            throw(DomainError(thickness,
                "air and bottom half-spaces must be infinite; internal layers positive and finite"))
        media(formula) === Val(:homogeneous) && length(thickness) != 2 &&
            throw(DimensionMismatch("homogeneous equations have no internal soil interfaces"))
    end
    return formula
end

"Route an explicit external formulation tag to its `Val` dispatch method."
Formulation(backend::Symbol; kwargs...) = Formulation(Val(backend); kwargs...)

function _fem_formulation(
        insulation_admittance, semicon_admittance, earth_properties, temperature_dependence,
        options::NamedTuple
)
    methods = (
        insulation_admittance = InsulationAdmittance.Formula(insulation_admittance),
        semicon_admittance = SemiconAdmittance.Formula(semicon_admittance),
        earth_properties = earth_properties === nothing ? nothing :
                           Earth.FrequencyDependent.Formula(earth_properties),
        temperature_dependence = temperature_dependence === nothing ? nothing :
                                 TemperatureDependent.Formula(temperature_dependence),
    )
    definitions = (; insulation_admittance, semicon_admittance, earth_properties,
        temperature_dependence)
    return LineCableModelsFEM(methods, formulation_options(LineCableModelsFEM, options),
        definitions)
end

"""
$(TYPEDSIGNATURES)

Construct the Gmsh/GetDP finite-element formulation. FEM owns its field equations
and selects four material laws. Each law and `options` accepts a
scalar or an explicit `Grid`/`Gridspace`; varying inputs return a
`Gridspace{LineCableModelsFEM}`.

# Keywords

- `insulation_admittance`: Insulation admittivity law; `:default` is lossless.
- `semicon_admittance`: Semicon admittivity law; `:default` is lossless.
- `earth_properties`: Soil frequency-dependent constitutive law; `:default` and
  `nothing` preserve the declared static soil. Equivalent-earth reductions are
  unsupported. Air uses its declared static properties.
- `temperature_dependence`: Cable-material resistivity law; `:default` selects
  the linear law and `nothing` retains reference resistivity. Operating
  temperature belongs to `LineParametersProblem`.
- `options=(;)`: Field model (`physics=:quasi_tem` or `:quasi_fw`) and bundle,
  Kron, and ideal-transposition reductions. Hyphenated strings and symbols
  are also accepted for `physics`. Julia parses `:quasi-fw` as subtraction;
  use `:quasi_fw` or `Symbol("quasi-fw")`.
  Pass execution controls to `compute(...; options=(...))`.
- `combine=:product`: Product or zip composition among varying inputs.

Analytical impedance/admittance kernel keywords are rejected. Supported enclosure
geometry is represented directly in the FEM domain.
"""
function Formulation(
        ::Val{:LineCableModelsFEM};
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
        earth_properties = formula(:default),
        temperature_dependence = formula(:default),
        options = (;),
        combine::Symbol = :product
)
    return parameterize(
        LineCableModelsFEM,
        _fem_formulation,
        (insulation_admittance, semicon_admittance, earth_properties,
            temperature_dependence, options);
        combine
    )
end

function LineCableModelsFEM(; kwargs...)
    return Formulation(Val(:LineCableModelsFEM); kwargs...)
end

function validate(binding::FormulaMethod, reduction::EquivalentHomogeneous.AbstractRule)
    throw(ArgumentError("$binding does not admit equivalent-earth reduction :$(formula_id(reduction))"))
end

"""Expose FEM constitutive/admittance selections, field model, and reductions."""
function Base.NamedTuple(value::LineCableModelsFEM)
    record = function (selected)
        selected === nothing && return nothing
        selected isa Symbol && return NamedTuple(formula(selected))
        selected isa NamedTuple && return map(record,selected)
        return NamedTuple(selected)
    end
    Record=NamedTuple{(:backend,:requested,:methods,:options),
        Tuple{Symbol,NamedTuple,NamedTuple,NamedTuple}}
    return Record((:fem,map(record,value.definitions),map(record,value.methods),value.options))
end
