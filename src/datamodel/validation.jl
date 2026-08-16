"""
    _with_kwdefaults(::Type{Part}, kwargs)

Merge the strict part constructor's declared keyword defaults with caller
keywords. Radial stacking is handled by group/spec materializers; this
validation layer receives numeric physical state.
"""
@inline function _with_kwdefaults(::Type{Part}, kwargs::NamedTuple) where {Part}
    defaults = Validation.keyword_defaults(Part)
    defaults === () && return kwargs
    values = defaults isa NamedTuple ? defaults :
        NamedTuple{Validation.keyword_fields(Part)}(defaults)
    return merge(values, kwargs)
end
