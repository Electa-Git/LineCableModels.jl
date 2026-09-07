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

Supertype for formulation values that select a complete numerical backend.
"""
abstract type AbstractFormulationBackend <: AbstractFormulation end

"""
$(TYPEDEF)

Supertype for backend-owned formulation option records.
"""
abstract type AbstractFormulationOptions end

"""
$(TYPEDEF)

Configure execution owned by the Gmsh/GetDP finite-element backend.

Physical geometry, materials, earth properties, frequencies, reductions, and
temperature remain properties of `LineParametersProblem` or the shared
line-parameter formulation options.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEMOptions <: AbstractFormulationOptions
    "Open the optional Gmsh graphical interface."
    ui::Bool
    "Generate field-map files for every frequency and basis terminal."
    plot_field_maps::Bool
    "Mesh selection policy, either `:reuse` or `:remesh`."
    mesh_policy::Symbol
    "Optional existing Gmsh mesh path."
    mesh_path::Union{Nothing, String}
    "Retain a successful run directory."
    keep_run_directory::Bool
    "Optional GetDP executable path."
    getdp_executable::Union{Nothing, String}
    "Gmsh message verbosity from 0 through 5."
    gmsh_verbosity::Int
    "GetDP message verbosity from 0 through 5."
    getdp_verbosity::Int
end

"""
$(TYPEDSIGNATURES)

Construct validated finite-element execution options.

# Keywords

- `ui=false`: Open the Gmsh graphical interface.
- `plot_field_maps=false`: Emit spatial field maps for every solve.
- `mesh_policy=:reuse`: Reuse a compatible mesh or generate one. Use
  `:remesh` to regenerate it unconditionally.
- `mesh_path=nothing`: Optional existing `.msh` file.
- `keep_run_directory=false`: Retain successful run artifacts.
- `getdp_executable=nothing`: Explicit GetDP executable path.
- `gmsh_verbosity=2`: Gmsh message verbosity from 0 through 5.
- `getdp_verbosity=2`: GetDP message verbosity from 0 through 5.

# Returns

- A validated [`LineCableModelsFEMOptions`](@ref) value.

# Errors

- `ArgumentError`: An option has an unsupported value or an empty path.
"""
function LineCableModelsFEMOptions(;
        ui::Bool = false,
        plot_field_maps::Bool = false,
        mesh_policy::Symbol = :reuse,
        mesh_path::Union{Nothing, AbstractString} = nothing,
        keep_run_directory::Bool = false,
        getdp_executable::Union{Nothing, AbstractString} = nothing,
        gmsh_verbosity::Integer = 2,
        getdp_verbosity::Integer = 2
)
    mesh_policy in (:reuse, :remesh) || throw(ArgumentError(
        "mesh_policy must be :reuse or :remesh; got $(repr(mesh_policy))",
    ))
    normalized_mesh_path = mesh_path === nothing ? nothing : String(mesh_path)
    normalized_getdp = getdp_executable === nothing ? nothing :
                       String(getdp_executable)
    normalized_mesh_path === "" && throw(ArgumentError("mesh_path cannot be empty"))
    normalized_getdp === "" && throw(ArgumentError(
        "getdp_executable cannot be empty",
    ))
    gmsh_verbosity in 0:5 || throw(ArgumentError(
        "gmsh_verbosity must be an integer from 0 through 5",
    ))
    getdp_verbosity in 0:5 || throw(ArgumentError(
        "getdp_verbosity must be an integer from 0 through 5",
    ))
    return LineCableModelsFEMOptions(
        ui,
        plot_field_maps,
        mesh_policy,
        normalized_mesh_path,
        keep_run_directory,
        normalized_getdp,
        Int(gmsh_verbosity),
        Int(getdp_verbosity)
    )
end

"""
$(TYPEDEF)

Select the Julia-native Gmsh/GetDP quasi-TEM finite-element backend.

`options` stores the shared LineCableModels reduction and temperature policy.
`execution` stores only Gmsh/GetDP execution controls.

$(TYPEDFIELDS)
"""
struct LineCableModelsFEM{M <: NamedTuple, O <: NamedTuple, D <: NamedTuple} <: AbstractFormulationBackend
    "Shared scientific formula selections, independent of FEM execution controls."
    methods::M
    "Shared line-parameter formulation options."
    options::O
    "Requested formula definitions retained for provenance."
    definitions::D
    "Finite-element execution options."
    execution::LineCableModelsFEMOptions
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

@required EarthImpedanceFormulation begin
    validate(::EarthImpedanceFormulation, ::EarthPair)
    validate(::EarthImpedanceFormulation, ::Integer)
end

@required EarthAdmittanceFormulation begin
    validate(::EarthAdmittanceFormulation, ::EarthPair)
    validate(::EarthAdmittanceFormulation, ::Integer)
end

"Return whether an earth formulation consumes homogeneous or stratified media."
media(::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation}) = Val(:homogeneous)

"""
$(TYPEDSIGNATURES)

Check an earth formula against the physical earth model before numerical
evaluation. Homogeneous formulas may consume an EHEM reduction downstream;
stratified formulas consume the physical horizontal layers directly.

# Arguments

- `formula`: Resolved earth-impedance or earth-admittance formula.
- `earth`: Validated static earth model, including the air layer.

# Returns

- The same `formula`.

# Errors

- Throws `ArgumentError` for vertical interfaces with a horizontal stratified
  formula, or a native exception for a formula-incompatible layer inventory.
"""
function validate(
        formula::Union{EarthImpedanceFormulation, EarthAdmittanceFormulation},
        earth::EarthModel
)
    validate(earth)
    earth.vertical_layers && media(formula) === Val(:stratified) &&
        throw(ArgumentError(
            "stratified earth-return formulas require horizontal earth interfaces; " *
            "the selected EarthModel has vertical interfaces"))
    validate(formula, length(earth.layers))
    return formula
end

"Route an explicit external formulation tag to its `Val` dispatch method."
Formulation(backend::Symbol; kwargs...) = Formulation(Val(backend); kwargs...)

function _fem_execution_options(options::LineCableModelsFEMOptions)
    return options
end

function _fem_execution_options(options::NamedTuple)
    return LineCableModelsFEMOptions(; options...)
end

function _fem_formulation(
        internal_impedance, insulation_impedance, earth_impedance,
        insulation_admittance, semicon_admittance, earth_admittance,
        earth_properties, equivalent_earth, pipe_impedance,
        options::NamedTuple,
        fem_options::Union{NamedTuple, LineCableModelsFEMOptions}
)
    physical = _line_formulation(internal_impedance, insulation_impedance,
        earth_impedance, insulation_admittance, semicon_admittance, earth_admittance,
        earth_properties, equivalent_earth, pipe_impedance, options)
    return LineCableModelsFEM(physical.methods, physical.options, physical.definitions,
        _fem_execution_options(fem_options))
end

"""
$(TYPEDSIGNATURES)

Construct the Gmsh/GetDP finite-element formulation through the code-first
`Formulation` grammar. Shared formula slots, `options`, and `fem_options` accept scalar values or
explicit [`Grid`](@ref LineCableModels.ParametricBuilder.Grid)/
[`Gridspace`](@ref LineCableModels.ParametricBuilder.Gridspace) sources; a
varying call returns a `Gridspace{LineCableModelsFEM}`.

# Keywords

- Formula slots have the same names and identifiers as [`Formulation`](@ref).
  Dielectric `:default` explicitly selects lossless admittivity. FEM retains its
  field equations instead of evaluating analytical impedance kernels; those
  selections are recorded with their backend treatment.
- `options=(;)`: Shared bundle, Kron, ideal-transposition, and temperature
  options accepted by the line-parameter engine.
- `fem_options=(;)`: A [`LineCableModelsFEMOptions`](@ref) value or the
  equivalent named tuple.
- `combine=:product`: `:product` or `:zip` composition between varying option
  bundles.
"""
function Formulation(
        ::Val{:LineCableModelsFEM};
        internal_impedance = formula(:default),
        insulation_impedance = formula(:default),
        earth_impedance = formula(:default),
        insulation_admittance = formula(:default),
        semicon_admittance = formula(:default),
        earth_admittance = formula(:default),
        earth_properties = formula(:default),
        equivalent_earth = formula(:default),
        pipe_impedance = formula(:default),
        options = (;),
        fem_options = (;),
        combine::Symbol = :product
)
    return parameterize(
        LineCableModelsFEM,
        _fem_formulation,
        (internal_impedance, insulation_impedance, earth_impedance,
            insulation_admittance, semicon_admittance, earth_admittance,
            earth_properties, equivalent_earth, pipe_impedance, options, fem_options);
        combine
    )
end

function LineCableModelsFEM(; kwargs...)
    return Formulation(Val(:LineCableModelsFEM); kwargs...)
end
