"""
$(TYPEDEF)

Select the Julia-native Gmsh/GetDP finite-element backend.

`options` stores the field-model choice, prescribed Γ and matrix reductions.
`options=(physics=:helmholtz,)` is the supported field model and the default.
Prescribe the longitudinal propagation constant with `options=(physics=:helmholtz, Γ=value)`.
The default is zero. A finite scalar applies at every frequency, and a vector
contains one value per frequency in problem order.
Execution controls belong to `compute(...; options=(...))` and are validated
by [`computation_options`](@ref) for `LineCableModelsFEM`.

The axial current excitation is one ampere. The model retains conduction
and displacement through ``κ=σ+jωε`` \\[S/m\\], where ``ω`` is angular
frequency \\[rad/s\\], ``σ`` conductivity \\[S/m\\], and ``ε`` permittivity \\[F/m\\].

The Helmholtz model uses phasors ``e^{jωt-Γz}`` and the exact substitution
``A_z=a``, ``A_t=Γb``, ``φ=Γv``. Here ``Γ`` is prescribed \\[1/m\\],
``a`` has units \\[T m\\], ``b`` \\[T m²\\], and ``v`` \\[V m\\].
In the media outside the equipotential metal terminals the equations are

```math
-\\nabla_t\\!\\cdot[\\nu(\\nabla_t a+\\Gamma^2 b)]
+j\\omega\\kappa a-\\Gamma^2\\kappa v=0,
\\qquad
C^*(\\nu Cb)+(j\\omega\\kappa-\\Gamma^2\\nu)b
+\\kappa\\nabla_t v-\\nu\\nabla_t a=0,
\\qquad
-\\nabla_t\\!\\cdot[\\kappa(j\\omega b+\\nabla_t v)]
+\\kappa(j\\omega a-\\Gamma^2v)=0.
```

Here ``\\nu=1/\\mu`` \\[m/H\\], ``μ`` is permeability \\[H/m\\],
``Cb=∂_x b_y-∂_y b_x`` and ``C^*h=(∂_y h,-∂_x h)``.
The finite-conductivity axial ``a/u`` block supplies both the series voltage
drop and the normalized transverse-current source. There is no independent
electric excitation. A tree gauge removes the gradient freedom of ``b``.
Its tangential trace vanishes on every terminal contour and the entire outer
boundary. ``v=0`` on both the air-side and earth-side outer boundaries.

``Γ^2`` is the complex square, not the squared modulus. These equations
retain all longitudinal terms. At exactly ``Γ=0`` they give the normalized
first-order limit. No small-Γ threshold or numerical division is used.
The exterior fields are ``E_t=-Γ(jωb+\\nabla_t v)`` and
``E_z=-jωa+Γ^2v``. Finite metal retains the axial model and the transverse
equipotential approximation. It solves for ``w_i=u_i-Γ^2v_i`` \\[V/m\\], giving
``E_z=-jωa-w_i`` inside terminal ``i``. This exact substitution cancels the
identical metal drive and potential basis contributions before assembly. It avoids
subtracting large terms to obtain the small driven current. The physical axial
drive is recovered as ``u_i=w_i+Γ^2v_i`` for output. At Γ=0, ``w_i=u_i``.

GetDP uses ``x`` horizontal, ``y`` vertical and ``z`` axial. The interface is
``y=0``. Each receiving row specifies the voltage reference for every source column:
air receivers use their local interface projection and buried receivers retain
earth infinity. The interface trace remains a solved field quantity.

For Helmholtz, paths ``ℓ_i`` run from the receiver reference to its
own conductor interior. The normalized inverse-admittance coefficient is

```math
P_{ij}=\\frac{\\int_{\\ell_i}(\\nabla v+j\\omega b)\\cdot d\\ell}{I_j},
\\qquad Y=P^{-1},\\qquad P_e=j\\omega P.
```

``I_j`` is the imposed axial current \\[A\\], ``P`` has units \\[Ω m\\],
``Y`` \\[S/m\\], and the analytical potential coefficient ``P_e`` \\[m/F\\].
The model forms complex primitive ``P`` before reductions and matrix inversion.
No additional ``jω`` factor is applied to ``Y``. Each terminal has one vertical
floating measurement curve ending inside its own terminal metal region.
The native horizontal shift is `min(1e-5,0.01*d)` \\[m\\], where `d` is the
terminal's smallest metal dimension \\[m\\]. Native GetDP integrates the stored
complex gauge-invariant vector ``\\nabla v+j\\omega b`` on media and PML.
There are no scalar endpoint terms or additional stretch factors.
Floating curves have their own nodes and do not enter the field support or
gauge tree. Stored transverse fields exclude metal interiors. Physical line sizing
uses one quarter of the local exterior native target. buried PML segments
use four times the bottom interval count with its progression. No division by ``Γ`` is used in extraction.
When comparing with the analytical `:unified` model, its mean-field receiver
(single line source) requires ``|κ_m r_p| ≲ 0.1``. Thick receiving
conductors near the interface are outside that declared comparison scope.
Exact air cutoff is singular, and a prescribed Γ can be non-passive.
The scalar trace ``v_i/I_j`` alone is gauge dependent and is saved separately
as `Pscalar.tsv` diagnostics. The longitudinal drive coefficient is
``K_{ij}=-U_i/I_j`` \\[Ω/m\\], where ``U_i`` is the axial drive \\[V/m\\].
The measured series coefficient is ``Z_{ij}=K_{ij}+Γ^2P_{ij}`` \\[Ω/m\\].
For a PEC contour, ``K_{ij}=(jωA_i-Γ^2v_i)/I_j``.

This is a two-dimensional Maxwell potential formulation with prescribed Γ,
transverse induction and displacement. It does not constitute the Darwin approximation.
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
    "Field model with prescribed longitudinal propagation constant and matrix reductions."
    options::O
    "Requested formula definitions."
    definitions::D
end

"""Identify the Gmsh/GetDP finite-element backend."""
description(::Type{<:LineCableModelsFEM}; compact::Bool = false) = "FEM"
function description(::LineCableModelsFEM; compact::Bool = false)
    description(LineCableModelsFEM; compact)
end
formula_id(::Type{<:LineCableModelsFEM}) = :fem
formula_id(::LineCableModelsFEM) = :fem
ImportExport.deserialize_value(::Val{:formulation}, ::Union{Val{:fem},Val{:LineCableModelsFEM}}) = LineCableModelsFEM
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

Report a failure in finite-element model preparation, meshing, solving, or result validation.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEMError <: Exception
    "Failure category."
    category::Symbol
    "Stable identifier of the object that records the failure."
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
scalar or an explicit `Grid`/`Gridspace`. varying inputs return a
`Gridspace{LineCableModelsFEM}`.

# Keywords

- `insulation_admittance`: Insulation admittivity law. `:default` routes to
  `:lossless`.
- `semicon_admittance`: Semicon admittivity law. `:default` routes to
  `:lossless`.
- `earth_properties`: Soil frequency-dependent constitutive law. `:default`
  routes to the explicit `:constant` pass-through, while `nothing` preserves
  the declared static soil. Equivalent-earth reductions are unsupported. Air
  uses its declared static properties.
- `temperature_dependence`: Cable-material resistivity law. `:default` selects
  the linear law and `nothing` retains reference resistivity. Operating
  temperature belongs to `LineParametersProblem`.
- `options=(;)`: Field model (`physics=:helmholtz`, the default), prescribed
  `Γ=0` [1/m] (scalar or frequency-aligned vector), and bundle,
  Kron, and ideal-transposition reductions. The string `"helmholtz"`
  is also accepted for `physics` and normalized to `:helmholtz`.
  Pass execution controls to `compute(...; options=(...))`.
- `combine=:product`: Product or zip composition among varying inputs.

Analytical impedance and admittance kernel keywords are rejected. Supported enclosure
geometry is represented directly in the FEM domain.
"""
function Formulation(
        ::Val{:LineCableModelsFEM};
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

function LineCableModelsFEM(; kwargs...)
    return Formulation(Val(:LineCableModelsFEM); kwargs...)
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

"""
Record the consumed FEM material laws and selected field assumptions.
"""
function formulation_record(formulation::LineCableModelsFEM)
    return merge((
        schema_version = 7,
        assumptions = (
            impedance = "Axial current-driven Maxwell equations; Z = -U/I + Gamma^2 P",
            admittance = "Coupled prescribed-Gamma Maxwell A_z/A_t/phi equations; axial current supplies normalized leakage; vertical path voltage includes A_t/Gamma; Y = inv(P)",
            earth = "Horizontal air and one semi-infinite soil; soil constitutive properties evaluated at each frequency",
            propagation = "Prescribed complex Gamma; exact A_t/Gamma and phi/Gamma variables with a regular zero limit; all Gamma^2 terms retained",
            semicon_domain = "Passive material region, without electrical terminal ownership",
            enclosure = "Supported enclosures are represented by their material and terminal domains"
        )
    ),NamedTuple(formulation))
end
