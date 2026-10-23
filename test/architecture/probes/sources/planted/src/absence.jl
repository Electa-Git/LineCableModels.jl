compared(x) = x === nothing || x !== nothing || nothing === x || ===(x, nothing)
equal(x) = x == nothing || x != missing || missing === x
tested(x) = isnothing(x) || Base.isnothing(x) || all(!isnothing, x) || any(isnothing, x)
missed(x) = ismissing(x) || something(x, 1) == 1
annotated(x::Nothing, ::Nothing; y::Nothing = nothing, z::Base.Nothing) = x
positional(x::Nothing = nothing) = x
absent(x::Missing; y::Union{Missing, Int} = missing) = x
anonymous = (x::Nothing) -> x
single = x::Nothing -> x
mapped(values) =
    map(values) do x::Nothing
        x
    end
struct Holder{T <: Union{Nothing, Int}}
    value::Union{Int, Nothing}
end
bounded(x::T) where {T <: Union{Real, Nothing}} = x
parametric(x::Union{Nothing, Int}, ::Type{Nothing}) = Quantity{Name, Core.Nothing}()
laundered(xs) = (i = first(xs); i === nothing)
mixed(xs) = (i = findfirst(isodd, xs); i = first(xs); i === nothing)
