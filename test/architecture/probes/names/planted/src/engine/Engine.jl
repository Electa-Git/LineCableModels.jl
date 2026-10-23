module Engine
import ..Commons
struct Formula <: Commons.AbstractFormulation end
struct Plan end
Commons.compute(::Formula) = 1
end
