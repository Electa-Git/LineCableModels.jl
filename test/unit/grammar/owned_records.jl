@testitem "Grammar / owned records / exact storage and value semantics" tags=[:unit] begin
    payload = (rtol=1e-8, nested=(enabled=true,), values=[1.0, 2.0], data=:field)
    @test (@inferred FormulationOptions(payload)).data === payload
    @test (@inferred ComputationOptions(payload)).data === payload
    @test (@inferred ComputationDetails(payload)).data === payload
    @test FormulationOptions !== ComputationOptions !== ComputationDetails
    @test LineCableModels.FormulaDefinition === LineCableModels.Grammar.FormulaDefinition
    @test LineCableModels.FormulaMethod === LineCableModels.Grammar.FormulaMethod

    for Record in (FormulationOptions, ComputationOptions, ComputationDetails)
        value = Record(payload)
        @test !(value isa NamedTuple)
        @test fieldnames(typeof(value)) == (:data,)
        @test fieldtype(typeof(value), :data) === typeof(payload)
        @test value.data === payload
        @test value.data.values === payload.values
        @test value.data.data === :field
        @test Record(; payload...).data === payload
        @test Record().data === NamedTuple()
        @test Record((other=3,)).data === (other=3,)
        @test_throws MethodError Record(Dict(:rtol=>1e-8))
        @test_throws MethodError Record((1, 2))
        @test_throws MethodError Record(3)
        @test_throws MethodError Record{NamedTuple}(payload)
        @test_throws MethodError Record(value)
        for Other in (FormulationOptions, ComputationOptions, ComputationDetails)
            @test_throws MethodError Record(Other(payload))
            Other === Record && continue
            @test value != Other(payload)
            @test !isequal(value, Other(payload))
        end
        @test value != payload
        @test !isequal(value, payload)
        for (a, b) in (((values=[1, 2],), (values=[1.0, 2.0],)),
                ((x=NaN,), (x=NaN,)), ((x=0.0,), (x=-0.0,)),
                ((x=missing,), (x=missing,)))
            @test isequal(Record(a) == Record(b), a == b)
            @test isequal(Record(a), Record(b)) == isequal(a, b)
            isequal(a, b) && @test hash(Record(a)) == hash(Record(b))
        end
        bounded = NamedTuple{(:shunt_model,), Tuple{NamedTuple}}(((iterations=4,),))
        @test fieldtype(typeof(Record(bounded).data), :shunt_model) === NamedTuple
        @test occursin(string(nameof(Record)), sprint(show, value))
    end

    retained = ComputationDetails(payload)
    payload.values[1] = 7.0
    @test retained.data.values[1] == 7.0
end

@testitem "Grammar / owned records / owner-specific option normalization" tags=[:unit] begin
    using LineCableModels.Grammar: formulation_options, computation_options
    resolved = @inferred formulation_options(LineParametersFormulation, FormulationOptions())
    @test resolved isa FormulationOptions
    @test resolved.data == (reduce_bundle=true, kron_reduction=true, ideal_transposition=true)
    execution = @inferred computation_options(LineCableModelsCoaxial, ComputationOptions())
    @test execution isa ComputationOptions
    @test execution.data.trace === Val(false)
    @test execution.data.output_basis === Val(:pul)
    @test_throws MethodError formulation_options(LineParametersFormulation, (;))
    @test_throws MethodError formulation_options(LineParametersFormulation, ComputationOptions())
    @test_throws MethodError computation_options(LineCableModelsCoaxial, (;))
    @test_throws MethodError computation_options(LineCableModelsCoaxial, FormulationOptions())
    @test_throws ArgumentError formulation_options(LineParametersFormulation, FormulationOptions(unknown=true))
    @test_throws ArgumentError computation_options(LineCableModelsCoaxial, ComputationOptions(trace="yes"))
    @test_throws ArgumentError computation_options(LineCableModelsCoaxial, ComputationOptions(ui=true))
    @test_throws TypeError formula(:default; options=ComputationDetails())
end

@testitem "Grammar / owned records / explicit transport and passive persistence" tags=[:unit] setup=[TestFixtures] begin
    using Serialization
    const IO = LineCableModels.ImportExport
    payload = (resolution=(enabled=true,), values=[1.0, 2.0])
    for Record in (FormulationOptions, ComputationOptions, ComputationDetails)
        record = Record(payload)
        @test_throws MethodError keys(record)
        @test_throws MethodError iterate(record)
        @test_throws MethodError NamedTuple(record)
        @test_throws MethodError merge(record, (;))
        for encoded in (IO.serialize_value(record), IO.serialize_value(record, Val(:scientific)))
            @test encoded["__type__"] == string(nameof(Record))
            restored = IO.deserialize_value(encoded)
            @test typeof(restored) === typeof(record)
            @test isequal(restored, record)
        end
        buffer = IOBuffer()
        serialize(buffer, record)
        seekstart(buffer)
        @test isequal(deserialize(buffer), record)
    end

    problem = CableConstantsProblem(TestFixtures.coaxial_design())
    formulation = CableConstantsFormulation()
    scalar = @inferred compute(problem, formulation; options=ComputationOptions())
    bundled = @inferred compute(problem, [formulation, formulation])
    @test isconcretetype(eltype(bundled))
    @test all(value -> typeof(value) === typeof(scalar), bundled)
    @test scalar == compute(problem, formulation; options=(;))
    @test details(scalar) isa ComputationDetails
    @test fieldtype(typeof(scalar), :details) === typeof(details(scalar))
    @test fieldtype(typeof(details(scalar).data), :shunt_model) === NamedTuple
    restored_scalar=IO.deserialize_value(IO.serialize_value(scalar))
    @test typeof(restored_scalar) === typeof(scalar)
    @test fieldtype(typeof(details(restored_scalar).data), :shunt_model) === NamedTuple
    @test_throws TypeError compute(problem, formulation; options=FormulationOptions())
    @test_throws MethodError ParametricProblem(problem, (;))

    for outer in (LinearError(formulation; options=(retain_details=true,)),
            MonteCarlo(formulation; trials=2, seed=13, retain_details=true, return_samples=true))
        result = compute(ParametricProblem(problem), outer)
        @test fieldtype(typeof(result), :details) === typeof(details(result))
        @test isconcretetype(eltype(result.values))
        encoded = IO.serialize_value(result)
        retained = IO.deserialize_value(encoded["details"])
        @test retained isa NamedTuple
        records = outer isa LinearError ? retained.points : first(retained.trials)
        @test all(record -> record isa NamedTuple, records)
        restored = IO.deserialize_value(encoded)
        @test isequal(details(restored), details(result))
        @test first(restored) == first(result)
        @test restored.formulation isa NamedTuple
        @test fieldtype(typeof(details(first(restored)).data), :shunt_model) === NamedTuple
        restored_records=outer isa LinearError ? details(restored).data.points :
            first(details(restored).data.trials)
        @test all(record -> fieldtype(typeof(record.data), :shunt_model) === NamedTuple,
            restored_records)
    end
end
