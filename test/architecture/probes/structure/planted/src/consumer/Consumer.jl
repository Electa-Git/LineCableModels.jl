module Consumer
import ..Early
export filter, count
struct Local end
function quantity end
function filter end
function count end
filter(::Local) = 1
Early.extend(::Int) = 1
Early.Item(::Int) = Early.Item()
include("../early/stray.jl")
include("../../ext/core.jl")
end
