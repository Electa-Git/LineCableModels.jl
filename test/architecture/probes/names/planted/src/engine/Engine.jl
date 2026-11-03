module Engine
import ..Commons
struct Formula <: Commons.AbstractFormulation
    Γ::Float64
    γ::Float64
end
struct Plan end
Commons.compute(::Formula; tolerance = 0.1) = 1
end
