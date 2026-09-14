# Fresh current-schema executable transport. Both processes define the same
# current builder; no substituted source or compatibility-era callable exists.
using LineCableModels, Measurements, JLD2, Test, Random
include(joinpath(@__DIR__,"../../../gauntlet/Gauntlet.jl"))
include(joinpath(@__DIR__,"../../support/scenarios.jl"))
using .Gauntlet, .CurrentScenarios
function transport_problem(scale, temperature)
    design=CurrentScenarios.coaxial_design(;scale)
    system=build(LineCableSystem,design,Pose2(.03,-1.0);
        connections=Dict(:core=>1,:sheath=>2))
    LineParametersProblem(system;temperature,earth_props=homogeneous(rho=100.0),
        frequencies=[10.0,1000.0])
end
mode,path=ARGS
if mode=="write"
    space=Gridspace{LineParametersProblem}(transport_problem,
        (Grid((1.0,1.1),AbsoluteError(.01)),Grid(20.0,AbsoluteError(1.0))))
    JLD2.jldsave(path;model_bytes=Gauntlet._execution_bytes(space))
elseif mode=="read"
    restored=Gauntlet._read_execution(path,"model")
    @testset "current executable inputs survive process boundary" begin
        @test length(restored)==2
        for (index,problem) in enumerate(restored)
            @test problem.system.designs[1].terminal_order==[:core,:sheath]
            @test problem.frequencies==[10.0,1000.0]
            @test nominal(problem.temperature)==20.0
            @test uncertainty(problem.temperature)==1.0
            @test nominal(outer_radius(problem.system.designs[1])) ≈ .012*(index==1 ? 1.0 : 1.1)
        end
        first_sample=rand(Xoshiro(2029),restored;distribution=:uniform)
        second_sample=rand(Xoshiro(2029),restored;distribution=:uniform)
        @test Gauntlet.numerical_input_sha256(first_sample)==Gauntlet.numerical_input_sha256(second_sample)
    end
else
    throw(ArgumentError("unknown transport operation: $mode"))
end
