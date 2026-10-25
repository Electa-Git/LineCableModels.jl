module FormulasProbeClean
module Commons
function defaults end
end
module Laws
import ..Commons: defaults
const FORMULAS = (:linear,)
struct Formula{ID} end
defaults(::Type{Formula{ID}}) where {ID} = throw(ArgumentError("unknown law :$ID"))
defaults(::Type{Formula{:linear}}) = (;)
rate(::Formula{:linear}, x) = x
end
module Kernels
const FORMULAS = ()
struct Formula{ID} end
end
end
