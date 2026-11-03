module FormulasProbePlanted
module Commons
function defaults end
function validate end
end
module Laws
import ..Commons: defaults, validate
const FORMULAS = (:linear, :idle, :default)
struct Formula{ID} end
defaults(::Type{Formula{ID}}) where {ID} = throw(ArgumentError("unknown law :$ID"))
defaults(::Type{Formula{:linear}}) = (;)
defaults(::Type{Formula{:idle}}) = (;)
defaults(::Type{Formula{:orphan}}) = (;)
defaults(::Type{Formula{:default}}) = (;)
validate(formula::Formula{:idle}) = formula
rate(::Formula{:linear}, x) = x
rate(::Formula{:orphan}, x) = 2x
(formula::Formula{:default})(x) = x
end
end
