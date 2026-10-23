module Alpha
import ..Commons
function uses(m)
    Commons.shared(1) + Commons.lonely(1) + Commons.undocumented(1) +
    Commons.untested(1) + Commons.ideal_transposition!(m)[1] +
    Commons.hidden_wrapper([1])[1] + Commons.extension_shared(1)
end
alone(::Commons.Lone) = 1
const EPS0 = 8.8541878128e-12
transposed!(m) = m
function permeability(T)
    μ0 = 4π * 1e-7
    return T(μ0)
end
_tiny(x) = x + 1
twice(x) = 2 * _tiny(x)
_named(x) = x + 2
uses_named(x) = 3 * _named(x)
Commons.extended(x::Int) = x
end
