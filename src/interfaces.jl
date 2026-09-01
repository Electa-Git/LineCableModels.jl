"Add one owned value to a mutable collection."
function add! end

"""
    build(Target, declarations...; kwargs...)

Construct one completed domain object from complete physical declarations.

`build(CableDesign, ...)` resolves a physical cable root and terminal state.
`build(LineCableSystem, ...)` places completed designs and resolves global
connections. When an argument is a `Grid` or `Gridspace`, or an admitted tuple
or vector contains one, the same surface returns a `Gridspace{Target}` whose
callable invokes scalar `build` after selecting and reconstructing one complete
point.

# Arguments

- `Target`: Completed domain type to construct.
- `declarations`: Complete physical declarations owned by `Target`.

# Returns

One completed `Target`, or `Gridspace{Target}` for explicit finite inputs.
"""
function build end

"""
    Gridpoint{Target}(build, args)

Store one selected but unresolved argument tuple from a finite parameter space.
`Target` preserves the semantic object or problem family that the point will
materialise, allowing computation dispatch to consume a scalar point without
first discarding its target identity.
"""
struct Gridpoint{Target, F, A <: Tuple}
    "Scalar constructor or lowering function selected by the finite space."
    build::F
    "Selected argument tuple, with any nested target-bearing points retained."
    args::A
end

function Gridpoint{Target}(build, args::A) where {Target, A <: Tuple}
    return Gridpoint{Target, typeof(build), A}(build, args)
end

# Shared scalar-versus-parametric construction boundary. Data-owning modules
# call `_construction` without depending on Gridspace. ParametricBuilder marks
# its finite sources and owns the lazy branch after it is loaded.
_construction_axis(::Any) = false

function _finite_construction end

function _construction(
        ::Type{Target},
        caller,
        values::Tuple;
        combine::Symbol = :product
) where {Target}
    combine in (:product, :zip) || throw(ArgumentError(
        "combine must be :product or :zip"
    ))
    any(_construction_axis, values) || return caller(values...)
    return _finite_construction(Target, caller, values; combine)
end

"""
    homogenize(design; new_id="", dielectric_frequency=50)

Build a homogeneous cable design that preserves each radial assembly member
and matches the effective conductor and dielectric properties of `design`.

The physical source design remains unchanged. The reduction uses only scalar
series and parallel circuit calculations; it does not calculate line-parameter
matrices, mutual coupling, or earth return.

# Arguments

- `design`: Completed physical cable design.

# Keywords

- `new_id`: Identifier for the returned design. An empty value appends
  `"_equivalent"` to the source identifier.
- `dielectric_frequency`: Reference frequency used to match lossy dielectric
  admittance \\[Hz\\]. Default: `50`.

# Returns

- A completed homogeneous `CableDesign`.
"""
function homogenize end

"Evaluate the constitutive relation selected for a material value."
function constitutive end

"Return the physical storage basis of a result."
function basis end

function R end
function L end
function C end
function resistance end
function inductance end
function capacitance end
