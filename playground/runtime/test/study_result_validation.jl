# Run in the actual line-parameter project, never inside a UI or coordinator.
# Match worker/test/scientific.jl's direct engine reference. Cross-machine
# cancellation residuals near zero require an absolute roundoff floor as well.
using Test, LineCableModelsLineParameters
const Adapter = LineCableModelsLineParameters
Adapter.load_profile!(Adapter.LineParameterProfile())
const Engine = Adapter.LineCableModels
const ExecutionCore = Adapter.LineCableModelsExecutionCore
length(ARGS)==1 || error("Expected the bounded browser result file")
filesize(ARGS[1]) <= 4*1024^2 || error("Scientific evidence exceeds its bound")
data = ExecutionCore.JSON3.read(read(ARGS[1],String),Dict{String,Any})

function matches_coefficient(actual, expected, scale)
    # One Float64 spacing at the frequency slice's coefficient scale; no fixed
    # tolerance with physical units and no relative tolerance against a zero.
    isapprox(actual, expected; rtol=sqrt(eps(Float64)), atol=eps(Float64)*scale)
end

@testset "cross-machine parity tolerance stays bounded" begin
    for scale in (1.0e-11, 1.0e-4)
        @test matches_coefficient(eps(Float64)*scale/2, 0.0, scale)
        @test !matches_coefficient(1.0e-8*scale, 0.0, scale)
        @test !matches_coefficient(scale*(1+1.0e-6), scale, scale)
    end
    @test !matches_coefficient(NaN, 0.0, 1.0)
    @test !matches_coefficient(eps(Float64), 0.0, 0.0)
end

@testset "both live consumers match the direct line-parameter engine" begin
    parameters = Adapter.validate_line_parameters(data["deck"]["parameters"]["parameters"])
    design = Adapter.build_coaxial_design(parameters)
    half = parameters["separation_m"]/2
    system = Engine.build(Engine.LineCableSystem,[design,design],
        [Engine.Pose2(-half,-parameters["depth_m"],0.0),Engine.Pose2(half,-parameters["depth_m"],0.0)];
        connections=[Dict(:core=>1,:sheath=>0),Dict(:core=>2,:sheath=>0)],
        system_id="registered-consumer-parity",line_length=parameters["line_length_m"])
    reference = Engine.compute(Engine.LineParametersProblem(system;
        temperature=parameters["temperature_celsius"],
        earth_props=Engine.Earth(rho=parameters["earth_resistivity_ohm_m"]),
        frequencies=parameters["frequencies_hz"]),Engine.Formulation())
    for consumer in ("deck","workbench")
        entry = data[consumer]["parameters"]
        @test entry["parameters"] == parameters
        @test entry["receipt"]["input_hash"] == Adapter.input_hash("line.frequency_scan",parameters)
        @test entry["value"]["frequencies_hz"] == parameters["frequencies_hz"]
        for (field,matrix) in (("series_impedance_ohm_per_m",Engine.Z(reference)),
                ("shunt_admittance_s_per_m",Engine.Y(reference)))
            for index in CartesianIndices(matrix)
                i,j,k = Tuple(index)
                value = entry["value"][field][i][j][k]
                scale = maximum(abs, view(matrix,:,:,k))
                @test matches_coefficient(ComplexF64(value["real"],value["imag"]), matrix[index], scale)
            end
        end
    end
    @test !any(id.name in ("Bonito","LineCableModelsPlayground","NATS","PowerImpedance")
        for id in keys(Base.loaded_modules))
end
