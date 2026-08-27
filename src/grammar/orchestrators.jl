"""
$(TYPEDSIGNATURES)

Return the abstract owner declared for a public action by [`@orchestrator`](@ref).

Missing declarations throw `MethodError`.
"""
function orchestrator_root end

"""
$(TYPEDSIGNATURES)

Return the sole root method declared for a public action and abstract owner by
[`@orchestrator`](@ref). The returned `Method` records its signature, module,
source file, and source line.

Missing declarations throw `MethodError`.
"""
function orchestrator_method end

function _declared_orchestrator_method(
        action::Function,
        owner::Module,
        file::Symbol,
        line::Integer
)
    candidates = Method[method
                        for method in methods(action)
                        if method.module === owner && method.file === file &&
                               method.line == line]
    length(candidates) == 1 || error(
        "expected one orchestrator method for $action at $file:$line; " *
        "found $(length(candidates))",
    )
    return only(candidates)
end

"""
    @orchestrator AbstractRoot function public_action(...)
        ...
    end

Define `public_action` unchanged and declare its method as the sole public
action for `AbstractRoot`. The declaration is available through
[`orchestrator_root`](@ref) and [`orchestrator_method`](@ref).

The macro accepts one ordinary long-form function definition. Hook
requirements remain owned by the definition hierarchy and may be declared
with `RequiredInterfaces.@required`.
"""
macro orchestrator(root, definition)
    definition isa Expr && definition.head === :function || throw(ArgumentError(
        "@orchestrator requires a long-form function definition",
    ))
    signature = definition.args[1]
    while signature isa Expr && signature.head === :where
        signature = signature.args[1]
    end
    signature isa Expr && signature.head === :call || throw(ArgumentError(
        "@orchestrator requires a named function definition",
    ))
    action = first(signature.args)
    action isa Symbol || throw(ArgumentError(
        "@orchestrator requires an unqualified function name",
    ))

    root_metadata = GlobalRef(@__MODULE__, :orchestrator_root)
    method_metadata = GlobalRef(@__MODULE__, :orchestrator_method)
    find_method = GlobalRef(@__MODULE__, :_declared_orchestrator_method)
    return quote
        #! explicit-imports: off
        # Julia's documentation protocol requires this non-public Base macro
        # when a macro emits the documented definition.
        Base.@__doc__ $(esc(definition))
        #! explicit-imports: on
        $root_metadata(
            ::typeof($(esc(action))),
            ::Type{$(esc(root))}
        ) = $(esc(root))
        $method_metadata(
            ::typeof($(esc(action))),
            ::Type{$(esc(root))}
        ) = $find_method(
            $(esc(action)),
            $(QuoteNode(__module__)),
            $(QuoteNode(__source__.file)),
            $(__source__.line)
        )
    end
end
