include("spectrum.jl")
include("fixtures.jl")

assessment=TOML.parsefile(joinpath(ROOT,"warmed-assessment.toml"))
haskey(assessment,"combined") && assessment["combined"]["relative_saving"]>0 ||
    error("Combined change has not established a warmed performance benefit")
spectrum!("combined")
if TOML.parsefile(joinpath(ROOT,"spectrum-combined.toml"))["passed"]
    fixtures!("combined")
else
    say("STOP combined change: spectrum preservation failed; no feature promotion")
end
