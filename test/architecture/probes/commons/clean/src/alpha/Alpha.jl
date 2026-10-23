module Alpha
import ..Commons
function capacitance(m, T)
    ε0 = Commons.vacuum_permittivity(T)
    tolerance = 1e-7
    return Commons.shared(ε0) * Commons.ideal_transposition!(m)[1] * tolerance
end
phase(x) = 2π * x
front_use(values) = Commons.front(values; factor = 2)
"The vacuum permittivity is 8.8541878128e-12 F/m."
documented(x) = x
_tested(x) = x + 1
calls_tested(x) = _tested(x) * 2
function _long(x)
    a = x + 1
    b = a * 2
    c = b - 3
    return c
end
calls_long(x) = _long(x) + 1
_reorder(x, y) = Commons.shared(y, x)
first_reorder(x) = _reorder(x, 1)
second_reorder(x) = _reorder(1, x)
struct Thing <: Commons.AbstractThing end
function Commons.accumulate_into!(values::Vector{Int})
    total = zero(eltype(values))
    for index in eachindex(values)
        total += values[index]
        values[index] = total
    end
    return values
end
end
