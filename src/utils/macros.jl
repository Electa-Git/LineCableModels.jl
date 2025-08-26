export @parameterize, @measurify

using MacroTools

macro parameterize(container_expr, union_expr)
    # Evaluate the Union type from the provided expression 
    local union_type
    try
        # Core.eval gets the *value* of the symbol passed in (e.g., the actual Union type)
        union_type = Core.eval(__module__, union_expr)
    catch e
        error("Expression `$union_expr` could not be evaluated. Make sure it's a defined const or type.")
    end

    # Sanity check
    if !(union_type isa Union)
        error("Second argument must be a Union type. Got a `$(typeof(union_type))` instead.")
    end

    # Base.uniontypes gets the component types, e.g., (Float64, Measurement{Float64})
    component_types = Base.uniontypes(union_type)

    # Define a recursive function to substitute the placeholder `_` 
    function substitute_placeholder(expr, replacement_type)
        # If the current part of the expression is the placeholder symbol,
        # we replace it with the target type (e.g., Float64).
        if expr == :_
            return replacement_type
            # If the current part is another expression (like `Array{_,3}`),
            # we need to recurse into its arguments to find the placeholder.
        elseif expr isa Expr
            # Rebuild the expression with the substituted arguments.
            new_args = [substitute_placeholder(arg, replacement_type) for arg in expr.args]
            return Expr(expr.head, new_args...)
            # Otherwise, it's a literal or symbol we don't need to change (e.g., `:Array` or `3`).
        else
            return expr
        end
    end

    # Build the list of new, concrete types 
    parameterized_types = [substitute_placeholder(container_expr, t) for t in component_types]

    # Wrap the new types in a single `Union{...}` expression and escape 
    final_expr = Expr(:curly, :Union, parameterized_types...)
    return esc(final_expr)
end

"""
    @measurify(function_definition)

Wraps a function definition. If any argument tied to a parametric type `T` is a
`Measurement`, this macro automatically promotes any other arguments of the same
parametric type `T` to `Measurement` with zero uncertainty. Other arguments
(e.g., `i::Int`) are forwarded without changes.
"""
macro measurify(def)
    # --- helpers at expand time ---
    function _contains_tvar(expr, typevars::Set{Symbol})
        found = false
        MacroTools.postwalk(expr) do x
            if x isa Symbol && x in typevars
                found = true
            end
            x
        end
        found
    end
    function _relax_container_type(ty, typevars::Set{Symbol}, bounds::Dict{Symbol,Any})
        if ty isa Expr && ty.head == :curly
            head = ty.args[1]
            params = Any[_contains_tvar(p, typevars) ?
                         MacroTools.postwalk(p) do x
                (x isa Symbol && haskey(bounds, x)) ? bounds[x] : x
            end |> x -> Expr(:<:, x) :
                         p
                         for p in ty.args[2:end]]
            return Expr(:curly, head, params...)
        else
            return MacroTools.postwalk(ty) do x
                (x isa Symbol && haskey(bounds, x)) ? bounds[x] : x
            end
        end
    end
    _ARRAY_HEADS = Set{Symbol}(
        [:Vector, :Array, :AbstractVector, :AbstractArray,
        :Matrix, :AbstractMatrix, :UnitRange, :StepRange, :AbstractRange]
    )
    _is_array_annot(ty) =
        ty isa Expr && ty.head == :curly && ty.args[1] isa Symbol && (ty.args[1] in _ARRAY_HEADS)
    # --- end helpers ---

    # Normalize and split
    if def.head == :(=)
        def = MacroTools.longdef(def)
    elseif def.head != :function
        error("@measurify must wrap a function definition")
    end
    dict = MacroTools.splitdef(def)

    where_items = get(dict, :whereparams, [])
    typevars = Set{Symbol}()
    bounds = Dict{Symbol,Any}()
    for w in where_items
        tv, ub = w isa Expr && w.head == :(<:) ? (w.args[1], w.args[2]) : (w, :Any)
        push!(typevars, tv)
        bounds[tv] = ub
    end

    posargs = get(dict, :args, [])
    kwargs = get(dict, :kwargs, [])

    # Collect symbols by role
    scalar_syms = Symbol[]     # positionals/keywords whose annot is exactly ::T
    array_syms = Symbol[]     # args like Vector{T}, AbstractArray{T}, ranges
    anchor_sym = nothing      # first parametric NON-array arg that contains T (e.g., EarthModel{T})

    # Build relaxed wrapper signature
    wrapper_pos = Expr[]
    wrapper_kw = Expr[]

    # Positionals
    for arg in posargs
        nm, ty, default = MacroTools.splitarg(arg)
        if (ty isa Symbol) && haskey(bounds, ty)
            push!(scalar_syms, nm)
            push!(wrapper_pos, Expr(:(::), nm, bounds[ty]))
        elseif (ty !== nothing) && _contains_tvar(ty, typevars)
            if _is_array_annot(ty)
                push!(array_syms, nm)
            elseif anchor_sym === nothing
                anchor_sym = nm
            end
            push!(wrapper_pos, Expr(:(::), nm, _relax_container_type(ty, typevars, bounds)))
        else
            push!(wrapper_pos, arg)
        end
    end

    # Keywords
    for kw in kwargs
        nm, ty, default = MacroTools.splitarg(kw)
        if (ty isa Symbol) && haskey(bounds, ty)
            push!(scalar_syms, nm)
        elseif (ty !== nothing) && _contains_tvar(ty, typevars) && _is_array_annot(ty)
            push!(array_syms, nm)
        end
        push!(wrapper_kw, MacroTools.postwalk(kw) do x
            (x isa Symbol && haskey(bounds, x)) ? bounds[x] : x
        end)
    end

    tight = def
    wrapper_dict = deepcopy(dict)
    wrapper_dict[:args] = wrapper_pos
    wrapper_dict[:kwargs] = wrapper_kw

    # ---- wrapper body ----
    # 1) TargetType: scalars-only; if there is an anchor (e.g., model::EarthModel{T}),
    #    use a WIDENING rule seeded with zero(T).
    scalar_vals = Any[:($s) for s in scalar_syms]
    target_decl = anchor_sym === nothing ?
                  :(TargetType = _coerce_args_to_T($(scalar_vals...))) :
                  quote
        Tanchor = first(typeof($(anchor_sym)).parameters)
        TargetType = isempty(($(scalar_vals...,))) ? Tanchor :
                     _coerce_args_to_T(zero(Tanchor), $(scalar_vals...))
    end

    # 2) If anchored, CONVERT the anchor only if its T differs (avoid clone-on-noop)
    anchor_convert = anchor_sym === nothing ? nothing : quote
        Tcur = first(typeof($(anchor_sym)).parameters)
        if Tcur !== TargetType
            $(anchor_sym) = convert(EarthModel{TargetType}, $(anchor_sym))
        end
    end

    # 3) Coerce arrays and scalars
    array_casts = [:($a = _coerce_array_to_T($a, TargetType)) for a in array_syms]
    scalar_casts = [:($s = _coerce_scalar_to_T($s, TargetType)) for s in scalar_syms]

    # 4) Forward call
    arg_names = [MacroTools.splitarg(a)[1] for a in posargs]
    kw_forwards = [Expr(:kw, MacroTools.splitarg(kw)[1], MacroTools.splitarg(kw)[1]) for kw in kwargs]
    forward_call = Expr(:call, dict[:name])
    !isempty(kw_forwards) && push!(forward_call.args, Expr(:parameters, kw_forwards...))
    append!(forward_call.args, arg_names)

    wrapper_body = quote
        $target_decl
        $(anchor_convert === nothing ? nothing : anchor_convert)
        $(array_casts...)
        $(scalar_casts...)
        $forward_call
    end
    # ----------------------

    # Drop where if no raw typevars remain in wrapper signature
    needs_where = any(_contains_tvar(arg, typevars) for arg in [wrapper_pos..., wrapper_kw...])
    if !needs_where
        delete!(wrapper_dict, :whereparams)
    end

    wrapper_dict[:body] = wrapper_body
    loose = MacroTools.combinedef(wrapper_dict)

    # Keep your doc footer
    return esc(Expr(:block, :(Base.@__doc__ $tight), loose))
end

"""
$(TYPEDSIGNATURES)

Automatically exports public functions, types, and modules from a module. This is meant for temporary development chores and should never be used in production code.

# Arguments

- None.

# Returns

- An `export` expression containing all public symbols that should be exported.

# Notes

This macro scans the current module for all defined symbols and automatically generates an `export` statement for public functions, types, and submodules, excluding built-in and private names. Private names are considered those starting with an underscore ('_'), as per standard Julia conventions.
	
# Examples

```julia
@autoexport
```
"""
macro autoexport()
    mod = __module__

    # Get all names defined in the module, including unexported ones
    all_names = names(mod; all=true)

    # List of names to explicitly exclude
    excluded_names = Set([:eval, :include, :using, :import, :export, :require])

    # Filter out private names (starting with '_'), module name, built-in functions, and auto-generated method symbols
    public_names = Symbol[]
    for name in all_names
        str_name = string(name)

        startswith(str_name, "@_") && continue  # Skip private macros
        startswith(str_name, "_") && continue  # Skip private names
        name === nameof(mod) && continue  # Skip the module's own name
        name in excluded_names && continue  # Skip built-in functions
        startswith(str_name, "#") && continue  # Skip generated method symbols (e.g., #eval, #include)

        if isdefined(mod, name)
            val = getfield(mod, name)
            if val isa Function || val isa Type || val isa Module
                push!(public_names, name)
            end
        end
    end

    return esc(Expr(:export, public_names...))
end