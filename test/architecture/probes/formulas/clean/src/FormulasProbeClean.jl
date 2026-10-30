module FormulasProbeClean
module Commons
function defaults end
end
module Laws
import ..Commons: defaults
const FORMULAS = (:linear, :cubic)
struct Formula{ID} end
defaults(::Type{Formula{ID}}) where {ID} = throw(ArgumentError("unknown law :$ID"))
defaults(::Type{Formula{:linear}}) = (;)
defaults(::Type{Formula{:cubic}}) = (;)
rate(::Formula{:linear}, x) = x
(formula::Formula{:cubic})(::Val{:J}, x) = x^3
end
module Kernels
const FORMULAS = ()
struct Formula{ID} end
end
end
