"""
$(TYPEDSIGNATURES)

Determines the promoted numeric element type for convenience constructors of component `C`. The promotion is computed across the values of `coercive_fields(C)`, extracted from the normalized `NamedTuple` `ntv` produced by [`validate!`](@ref). This ensures all numeric fields that participate in calculations share a common element type (e.g., `Float64`, `Measurement{Float64}`).

# Arguments

- `::Type{C}`: Component type \\[dimensionless\\].
- `ntv`: Normalized `NamedTuple` returned by `validate!` \\[dimensionless\\].
- `_order::Tuple`: Ignored by this method; present for arity symmetry with `_coerced_args` \\[dimensionless\\].

# Returns

- The promoted numeric element type \\[dimensionless\\].

# Examples

```julia
Tp = $(FUNCTIONNAME)(Tubular, (radius_in=0.01, radius_ext=0.02, material_props=mat, temperature=20.0), ())
```
"""
@inline _promotion_T(::Type{C}, ntv, _order::Tuple) where {C} =
    resolve_T((getfield(ntv, k) for k in coercive_fields(C))...)

"""
$(TYPEDSIGNATURES)

Builds the positional argument tuple to feed the **typed core** constructor, coercing only the fields returned by `coercive_fields(C)` to type `Tp`. Non‑coercive fields (e.g., integer flags) are passed through unchanged. Field order is controlled by `order` (a tuple of symbols), typically `(required_fields(C)..., keyword_fields(C)...)`.

# Arguments

- `::Type{C}`: Component type \\[dimensionless\\].
- `ntv`: Normalized `NamedTuple` returned by `validate!` \\[dimensionless\\].
- `Tp`: Target element type for numeric coercion \\[dimensionless\\].
- `order::Tuple`: Field order used to assemble the positional tuple \\[dimensionless\\].

# Returns

- A `Tuple` of arguments in the requested order, with coercions applied where configured.

# Examples

```julia
args = $(FUNCTIONNAME)(Tubular, ntv, Float64, (:radius_in, :radius_ext, :material_props, :temperature))
```
"""
@inline _coerced_args(::Type{C}, ntv, Tp, order::Tuple) where {C} =
    tuple((
        let k = s, v = getfield(ntv, s)
            (s in coercive_fields(C)) ? coerce_to_T(v, Tp) : v
        end
        for s in order
    )...)

"""
$(TYPEDSIGNATURES)

Utility for the constructor macro to *materialize* input tuples from either:

- A tuple literal expression (e.g., `(:a, :b, :c)`), or
- A bound constant tuple name (e.g., `_REQ_TUBULAR`).

Used to keep macro call sites short while allowing both styles.

# Arguments

- `mod`: Module where constants are resolved \\[dimensionless\\].
- `x`: Expression or symbol representing a tuple \\[dimensionless\\].

# Returns

- A standard Julia `Tuple` (of symbols or defaults).

# Errors

- `ErrorException` if `x` is neither a tuple literal nor a bound constant name.

# Examples

```julia
syms = $(FUNCTIONNAME)(@__MODULE__, :( :a, :b ))
syms = $(FUNCTIONNAME)(@__MODULE__, :_REQ_TUBULAR)
```
"""
_ctor_materialize(mod, x) = x === :(()) ? () :
                            x isa Expr && x.head === :tuple ? x.args :
                            x isa Symbol ? Base.eval(mod, x) :
                            error("@_ctor: expected tuple literal or const tuple, got $(x)")

using MacroTools: postwalk
"""
$(TYPEDSIGNATURES)

Generates a weakly‑typed convenience constructor for a component `T`. The generated method:

1. Accepts exactly the positional fields listed in `REQ`.
2. Accepts keyword arguments listed in `OPT` with defaults `DEFS`.
3. Calls `validate!(T, ...)` forwarding **variables** (not defaults),
4. Computes the promotion type via `_promotion_T(T, ntv, order)`,
5. Coerces only `coercive_fields(T)` via `_coerced_args(T, ntv, Tp, order)`,
6. Delegates to the numeric core `T(...)` with the coerced positional tuple.

`REQ`, `OPT`, and `DEFS` can be provided as tuple literals or as names of bound constant tuples. `order` is implicitly `(REQ..., OPT...)`.

# Arguments

- `T`: Component type (bare name) \\[dimensionless\\].
- `REQ`: Tuple of required positional field names \\[dimensionless\\].
- `OPT`: Tuple of optional keyword field names \\[dimensionless\\]. Defaults to `()`.
- `DEFS`: Tuple of default values matching `OPT` \\[dimensionless\\]. Defaults to `()`.

# Returns

- A method definition for the weakly‑typed constructor.

# Examples

```julia
const _REQ_TUBULAR = (:radius_in, :radius_ext, :material_props)
const _OPT_TUBULAR = (:temperature,)
const _DEFS_TUBULAR = (T₀,)

@_ctor Tubular _REQ_TUBULAR _OPT_TUBULAR _DEFS_TUBULAR

# Expands roughly to:
# function Tubular(radius_in, radius_ext, material_props; temperature=T₀)
#   ntv = validate!(Tubular, radius_in, radius_ext, material_props; temperature=temperature)
#   Tp  = _promotion_T(Tubular, ntv, (:radius_in, :radius_ext, :material_props, :temperature))
#   args = _coerced_args(Tubular, ntv, Tp, (:radius_in, :radius_ext, :material_props, :temperature))
#   return Tubular(args...)
# end
```

# Notes

- Defaults supplied in `DEFS` are **escaped** into the method signature (evaluated at macro expansion time).
- Forwarding into `validate!` always uses *variables* (e.g., `temperature=temperature`), never literal defaults.
- The macro is hygiene‑aware; identifiers `validate!`, `_promotion_T`, `_coerced_args`, and the type name are properly escaped.

# Errors

- `ErrorException` if `length(OPT) != length(DEFS)`.
"""
macro _ctor(T, REQ, OPT=:(()), DEFS=:(()))
    mod = __module__
    req = Symbol.(_ctor_materialize(mod, REQ))
    opt = Symbol.(_ctor_materialize(mod, OPT))
    dfx = _ctor_materialize(mod, DEFS)
    length(opt) == length(dfx) || error("@_ctor: OPT and DEFS length mismatch")

    # A) signature defaults (escape defaults)
    sig_kws = [Expr(:kw, opt[i], esc(dfx[i])) for i in eachindex(opt)]
    #    forwarding kwargs (variables, not defaults)
    pass_kws = [Expr(:kw, s, s) for s in opt]

    # B) flat order tuple
    order_syms = (req..., opt...)
    order = Expr(:tuple, (QuoteNode.(order_syms))...)

    ex = isempty(sig_kws) ? quote
        function $(T)($(req...))
            ntv = validate!($(T), $(req...))
            Tp = _promotion_T($(T), ntv, $order)
            local __args__ = _coerced_args($(T), ntv, Tp, $order)
            return $(T)(__args__...)
        end
    end : quote
        function $(T)($(req...); $(sig_kws...))
            ntv = validate!($(T), $(req...); $(pass_kws...))   # C) pass vars
            Tp = _promotion_T($(T), ntv, $order)
            local __args__ = _coerced_args($(T), ntv, Tp, $order)
            return $(T)(__args__...)
        end
    end

    # hygiene stays as you had it
    free = Set{Symbol}([:validate!, :_promotion_T, :_coerced_args, T])
    ex2 = postwalk(ex) do node
        node isa Symbol && (node in free) ? esc(node) : node
    end
    return ex2
end