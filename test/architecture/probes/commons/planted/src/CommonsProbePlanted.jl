module CommonsProbePlanted
include("commons/Commons.jl")
using .Commons: welcome
export welcome
include("alpha/Alpha.jl")
include("beta/Beta.jl")
root_helper() = 1
const ROOT_LIMIT = 2
end
