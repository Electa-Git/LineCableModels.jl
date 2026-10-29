validate(x::Int) = nothing
validate(::Float64) = 1
function validate(x::String)
    if isempty(x)
        throw(ArgumentError("empty"))
    end
end
function validate(x::Vector{Int})
    isempty(x) && return
    return x
end
Owner.validate(x::Char) = normalize(x)
function validate(x::Vector, context)
    validate(context)
    return x
end
_check_input(x) = x
function ensure_ready end
Owner.require_kind(x) = x
__validate_shape(x) = x
