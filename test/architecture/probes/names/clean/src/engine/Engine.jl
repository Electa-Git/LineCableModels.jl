module Engine
import ..Commons
struct Formula <: Commons.AbstractFormulation
    Γ::Float64
    γ::Float64
end
Commons.compute(::Formula; tolerance = 0.1) = 2
Λ(formula::Formula, ::Val{:self}, κ, r, h) = κ * r * h
end
