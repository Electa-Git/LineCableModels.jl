module Late
import ..Early
import Base: sum
export sum
struct Local end
sum(::Local) = 0
Early.extend(::Local) = 2
Early.quantity(::Local) = 3
Early.Item(::Local) = Early.Item()
downward() = Early.extend(Early.Item())
include("helpers.jl")
end
