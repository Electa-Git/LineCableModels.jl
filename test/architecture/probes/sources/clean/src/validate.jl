validate(x::Int) = x
function Owner.validate(x::Real, context)
    x > 0 || throw(ArgumentError("nonpositive"))
    foreach(context) do y
        return y
    end
    return x
end
validate(x::Symbol) = isempty(string(x)) ? error("empty") : x
function validate(x::String)
    if isempty(x)
        throw(ArgumentError("empty"))
    elseif length(x) > 9
        return x
    else
        x
    end
end
function validate(x::Char)
    try
        isvalid(x) || throw(ArgumentError("invalid"))
        x
    catch
        rethrow()
    end
end
validate(x::Float32, ::Type{Int}) = x
function validate end
validate!(x) = nothing
validated(x) = nothing
checker(x) = x
check(x) = x
call() = _check_input(1)
