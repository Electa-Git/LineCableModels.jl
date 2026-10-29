module Beta
import ..Commons
function admittance(m, T)
    Commons.shared(1) * Commons.vacuum_permittivity(T) *
    Commons.ideal_transposition!(m)[1]
end
Commons.accumulate_into!(values::Vector{Float64}) = values
front_use(values) = Commons.front(values; factor = 3)
end
