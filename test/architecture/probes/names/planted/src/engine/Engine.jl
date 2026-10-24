module Engine
import ..Commons
struct Formula <: Commons.AbstractFormulation end
struct Plan end
Commons.compute(::Formula; tolerance = 0.1) = 1
end
