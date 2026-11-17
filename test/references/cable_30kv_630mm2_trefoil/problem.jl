# The problem of the example `examples/cable_30kv_630mm2_trefoil.jl`: three 18/30 kV cables
# of 630 mm² in a touching trefoil. Its frequencies are the 50 multiples of 50 Hz up to
# 2500 Hz. The example runs in a module of its own.
let example = Module()
    Base.include(example,
        joinpath(@__DIR__, "..", "..", "..", "examples", "cable_30kv_630mm2_trefoil.jl"))
    Base.invokelatest(getglobal, example, :problem)
end
