returned() = nothing
defaulted(x; fallback = nothing) = x
located(xs, x) = findfirst(==(x), xs) === nothing
function bound(xs, x)
    i = findfirst(==(x), xs)
    i === nothing && throw(ArgumentError("absent"))
    return i
end
function walked(itr)
    next = iterate(itr)
    while next !== nothing
        (value, state) = next
        next = iterate(itr, state)
    end
end
looped(itr) =
    while (next = iterate(itr)) !== nothing
    end
parsed(text) = something(tryparse(Int, text), 0)
matched(text) = isnothing(match(r"a", text))
compared(x) = x !== Some(1)
annotated()::Nothing = nothing
text = "x === nothing, isnothing(x), ::Nothing, Union{Nothing, Int}" # isnothing(x)
