selected(xs) = [x for x in xs if x isa Int]
generated(xs) = (x for x in xs if typeof(x) <: Real)
filtered(xs) = filter(x -> x isa Int, xs)
found(xs) = findfirst(x -> x isa Int, xs)
located(xs) = Iterators.filter(x -> x isa Int, xs)
blocked(xs) =
    filter(xs) do x
        x isa Int
    end
function guarded(f)
    try
        f()
    catch e
        e isa ArgumentError || rethrow()
    end
end
bounded(x::T) where {T <: Real} = x
where_bound(x::T) where {T <: Real} = x
struct Wrapped{T <: Real} <: Number
    value::T
end
abstract type Kind <: Number end
typed(::Type{<:Real}) = 1
