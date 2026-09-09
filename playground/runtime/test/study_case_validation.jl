# Run separately in each actual scientific profile environment. The UI package
# must not be imported to compare its passive input contract with worker validation.
using Test
include(joinpath(@__DIR__,"..","..","src","scientific","StudyCases.jl"))
using .StudyCases

if only(ARGS) == "line-parameters"
    using LineCableModelsLineParameters
    const Adapter = LineCableModelsLineParameters
    @testset "line view inputs match authoritative worker normalization" begin
        for kwargs in ((;),(separation_m=0.8, depth_m=1.5, frequency_points=3),
                (earth_resistivity_ohm_m=500, minimum_frequency_hz=50, maximum_frequency_hz=1000))
            value = inputs(LineParameters(); kwargs...)
            normalized = Adapter.validate_line_parameters(value)
            @test value == normalized
            @test Adapter.input_hash(operation(LineParameters()),value) == Adapter.input_hash(operation(LineParameters()),normalized)
        end
        @test preparation_inputs(LineParameters()) == Adapter.validate_line_parameters(preparation_inputs(LineParameters()))
    end
elseif only(ARGS) == "power-flow"
    using LineCableModelsPowerFlow
    const Adapter = LineCableModelsPowerFlow
    @testset "corridor view inputs match authoritative worker normalization" begin
        for kwargs in ((;),(ugc_share=0.75, corridor_length_m=80000, frequency_points=3),
                (length_error_percent=0, minimum_frequency_hz=50, maximum_frequency_hz=1000))
            value = inputs(CorridorImpedance(); kwargs...)
            normalized = Adapter.validate_impedance_evaluation(value)
            @test value == normalized
            @test Adapter.input_hash(operation(CorridorImpedance()),value) == Adapter.input_hash(operation(CorridorImpedance()),normalized)
        end
        prepared = preparation_inputs(CorridorImpedance())
        @test prepared == Adapter.validate_powerflow_spec(prepared)
    end
else
    error("expected installed scientific profile name")
end
@test !any(id.name in ("Bonito","LineCableModelsPlayground","NATS") for id in keys(Base.loaded_modules))
@test !any(id.name in ("LineCableModels","PowerImpedance") for id in keys(Base.loaded_modules))
