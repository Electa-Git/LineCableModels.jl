module Engine
import ..Commons
struct Formula <: Commons.AbstractFormulation end
Commons.compute(::Formula; tolerance = 0.1) = 2
end
