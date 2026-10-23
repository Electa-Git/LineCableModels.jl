module Commons
export shared, vacuum_permittivity, ideal_transposition!, greet, accumulate_into!
export front, kernel, scale_of, backend, AbstractThing
public declared
"Applied inside the anonymous function of `front`."
kernel(x) = x + 1
"Applied inside the keyword body of `front`."
scale_of(x) = 2x
"Used by two owners. Its anonymous function and keyword body use the other two."
function front(values; factor = 1)
    weight = scale_of(factor)
    return map(value -> kernel(value) * weight, values)
end
"Used by two owners."
shared(x) = x + 1
"User API that no owner uses."
greet(x) = x
"User API that the root declares public."
declared(x) = x
"Grammar whose default method a method of `Alpha` repeats."
function accumulate_into!(values::AbstractVector)
    total = zero(eltype(values))
    for index in eachindex(values)
        total += values[index]
        values[index] = total
    end
    return values
end
"Grammar that only the extension extends."
function backend end
"Grammar that one module subtypes."
abstract type AbstractThing end
include("consts.jl")
include("matrixops.jl")
end
