import ..Validation: is_radius_input

"""
$(TYPEDSIGNATURES)

Default policy for **inner** radius raw inputs: accept proxies that expose an outer radius. This permits stacking by hijacking `p.radius_ext` during parsing.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `::Val{:radius_in}`: Field tag for the inner radius \\[dimensionless\\].
- `p::AbstractCablePart`: Proxy object \\[dimensionless\\].

# Returns

- `Bool` indicating acceptance (`true` if `hasproperty(p, :radius_ext)`).

# Examples

```julia
Validation.is_radius_input(Tubular, Val(:radius_in), prev_layer)  # true if prev_layer has :radius_ext
```
"""
is_radius_input(::Type{T}, ::Val{:radius_in}, p::AbstractCablePart) where {T} = hasproperty(p, :radius_ext)

"""
$(TYPEDSIGNATURES)

Default policy for **outer** radius raw inputs (annular shells): reject `AbstractCablePart` proxies. Outer radius must be numeric or a `Thickness` wrapper to avoid creating zero‑thickness layers.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `::Val{:radius_ext}`: Field tag for the outer radius \\[dimensionless\\].
- `::AbstractCablePart`: Proxy object \\[dimensionless\\].

# Returns

- `false` always.

# Examples

```julia
Validation.is_radius_input(Tubular, Val(:radius_ext), prev_layer)  # false
```
"""
is_radius_input(::Type{T}, ::Val{:radius_ext}, ::AbstractCablePart) where {T} = false

"""
$(TYPEDSIGNATURES)

Default policy for **outer** radius raw inputs (annular shells): accept `Thickness` as a convenience wrapper. The thickness is expanded to an outer radius during parsing.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `::Val{:radius_ext}`: Field tag for the outer radius \\[dimensionless\\].
- `::Thickness`: Thickness wrapper \\[dimensionless\\].

# Returns

- `Bool` indicating acceptance (`true`).

# Examples

```julia
Validation.is_radius_input(Tubular, Val(:radius_ext), Thickness(1e-3))  # true
```
"""
is_radius_input(::Type{T}, ::Val{:radius_ext}, ::Thickness) where {T} = true

"""
$(TYPEDSIGNATURES)

Merge per-part keyword defaults declared via `Validation.keyword_defaults` with
user-provided kwargs and return a **NamedTuple** suitable for forwarding.

Defaults may be a `NamedTuple` or a `Tuple` zipped against `Validation.keyword_fields(::Type{C})`.
User keys always win.
"""
@inline function _with_kwdefaults(::Type{C}, kwargs::NamedTuple) where {C}
    defs = Validation.keyword_defaults(C)
    defs === () && return kwargs
    nt = defs isa NamedTuple ? defs :
         NamedTuple{Validation.keyword_fields(C)}(defs)
    return merge(nt, kwargs)
end