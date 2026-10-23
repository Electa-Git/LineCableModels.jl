module Engine
import ..Commons
struct Formula <: Commons.AbstractFormulation end
Commons.compute(::Formula) = 2
end
