# Separate the already-qualified pilot changes on the first failing cable.
# Reuse all reference meshes/results; stop each candidate at its first failure.
include("fixtures.jl")

cases=[("screen",.1),("screen",1e6),("tube",.1),("tube",1e6),
    ("three",.1),("three",1e6),("sector",.1),("sector",1e6)]
for variant in ("harmonic","physical3")
    fixtures!(variant;cases,stop_on_failure=true)
end
say("COMPLETE separate assembly qualification; production equations unchanged")
