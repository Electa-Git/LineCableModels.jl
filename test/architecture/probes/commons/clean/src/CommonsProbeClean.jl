module CommonsProbeClean
include("commons/Commons.jl")
using .Commons: greet, declared
export greet
public declared
include("alpha/Alpha.jl")
include("beta/Beta.jl")
end
