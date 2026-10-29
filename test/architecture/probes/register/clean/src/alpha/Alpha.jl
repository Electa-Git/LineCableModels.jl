module Alpha
import ..Commons
public Item, area
struct Item{T <: Real}
    x::T
    Item{T}(x::T) where {T <: Real} = new{T}(x)
end
Item(x::Real) = Item{typeof(x)}(x)
(item::Item)(y) = item.x * y
Base.length(::Item) = 1
Commons.measure(item::Item) = item.x
area(item::Item) = item.x^2
@inline scale(item::Item, k) = item.x * k
magnitude(x::BigFloat) = sqrt(x * x)
magnitude(x) = abs(x)
const summed = sum
function totals(items::Vector{<:Item})
    total = sum(items; init = 0) do item
        area(item)
    end
    sorted = sort(items; by = item -> item.x)
    first = (; value = item -> item.x)
    return total + length(sorted) + first.value(items[1]) +
           sum(map(item -> scale(item, 2), items))
end
function __init__() end
end
