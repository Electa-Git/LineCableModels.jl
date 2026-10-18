@testitem "Linear error / derivative methods / owned values and construction" tags=[:unit,:uq] begin
    using LineCableModels.UQ: LinearErrorMethod, MeasurementsPropagation, CentralDifferences
    inner = CableConstantsFormulation()
    default = LinearError(inner)
    @test default.derivative isa MeasurementsPropagation
    @test default.options.data === (retain_details=false,)
    @test LinearError(inner, MeasurementsPropagation()).derivative === default.derivative
    @test LinearError(inner;derivative=MeasurementsPropagation()).derivative === default.derivative
    @test CentralDifferences().relative_step === 0.01
    @test CentralDifferences(.005).relative_step === .005
    @test CentralDifferences(relative_step=.005) === CentralDifferences(.005)
    @test LinearErrorMethod(:measurements) === MeasurementsPropagation()
    @test LinearErrorMethod(:central,.005) === CentralDifferences(.005)
    @test LinearErrorMethod(CentralDifferences(.005)) === CentralDifferences(.005)
    @test sprint(show,MeasurementsPropagation()) == "MeasurementsPropagation()"
    @test occursin("relative_step",sprint(show,CentralDifferences(.005)))
    @test occursin("0.005",sprint(show,MIME"text/plain"(),CentralDifferences(.005)))
    for method in (MeasurementsPropagation(),CentralDifferences(.005))
        direct = LinearError(inner,method;options=(retain_details=true,))
        keyword = LinearError(inner;derivative=method,retain_details=true)
        @test direct.derivative === method
        @test keyword.derivative === method
        @test isconcretetype(fieldtype(typeof(direct),:derivative))
        @test direct.options.data === (retain_details=true,)
        @test NamedTuple(direct) == NamedTuple(keyword)
        @test_throws ArgumentError LinearError(inner;derivative=method,relative_step=.01)
        @test_throws ArgumentError LinearError(inner,method;options=(unused=true,))
        @test_throws ArgumentError LinearError(inner,method,ComputationOptions(retain_details=1))
    end
    selected = LinearError(inner;derivative=:central,relative_step=.005,retain_details=true)
    stored = LinearError(inner;options=(derivative=:central,relative_step=.005,retain_details=true))
    direct = LinearError(inner,CentralDifferences(.005);options=(retain_details=true,))
    @test typeof(selected) === typeof(direct)
    @test selected.derivative === direct.derivative
    @test stored.derivative === direct.derivative
    @test NamedTuple(selected) == NamedTuple(direct) == NamedTuple(stored)
    @test NamedTuple(selected).options === (retain_details=true,derivative=:central,relative_step=.005)
    @test NamedTuple(default).options === (retain_details=false,derivative=:measurements,relative_step=.01)
    for step in (0,-1,Inf,NaN,true,nothing,"0.01")
        @test_throws ArgumentError CentralDifferences(step)
        @test_throws ArgumentError LinearError(inner;relative_step=step)
    end
    @test LineCableModels.validate(-1.0,CentralDifferences()) === -1.0
    @test_throws ArgumentError LineCableModels.validate(-1.0,CentralDifferences)
    @test_throws ArgumentError LinearError(inner,CentralDifferences();options=(relative_step=.01,))
    @test_throws ArgumentError LinearErrorMethod(:unknown)
end
