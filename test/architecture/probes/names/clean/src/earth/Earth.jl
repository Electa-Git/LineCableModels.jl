module Earth
import ..Commons
struct EarthModel
    layers::Vector{Float64}
    mu_r::Float64
end
abstract type AbstractLayer end
struct Pose2
    angle::Float64
end
struct LineSystem
    count::Int
end
struct Formula <: Commons.AbstractFormulation end
line_length(::LineSystem) = 1.0
figurelegend(x) = x
is_open(x) = true
stranded(x; compact = false) = x
homogenize(x) = x
Commons.compute(::Formula) = 1
Base.length(::LineSystem) = 2
end
