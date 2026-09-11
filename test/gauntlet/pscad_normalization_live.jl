# Explicit station check; this script is never discovered as a package test.
# Run from the repository root with --project=gauntlet.
using LineCableModels, JLD2, Test, EzXML

const P = LineCableModels.PSCAD
const E = LineCableModels.Engine
const OUTPUT_PARENT = joinpath(pkgdir(LineCableModels), ".linecablemodels",
    "gauntlet", "verification", "pscad-normalization")
mkpath(OUTPUT_PARENT)
const OUTPUT = mktempdir(OUTPUT_PARENT; prefix = "run-", cleanup = false)
cp(@__FILE__, joinpath(OUTPUT, basename(@__FILE__)))
length(ARGS) == 1 || error("supply an explicit station configuration file (see gauntlet/local.example)")
const OPTIONS = Base.include(Main,abspath(only(ARGS))).reference_options


function problem(; temperature = 60, rho = 1.72e-8, semicon = false)
    copper = Material(:conductor, rho, 1, 1, 20, 0.004)
    dielectric = Material(:insulator, 1e8, 2.3; tan_delta = 0.03)
    metal = solid(copper, Disk(0.004))
    outer = insulation(dielectric; t = 0.002)
    layers = semicon ? (metal,
        screen(Material(:semicon, 1e4, 40; tan_delta = 0.02); t = 0.0005), outer) :
        (metal, outer)
    design = build(CableDesign, "native-normalization", terminal(:core, layers...))
    system = build(LineCableSystem, fill(design, 4),
        [Pose2(0, 2), Pose2(1, 3), Pose2(2, -1), Pose2(3, -2)];
        connections = [Dict(:core => index) for index in 1:4],
        system_id = "native-normalization")
    return LineParametersProblem(system; temperature,
        earth_props = homogeneous(rho = 100.0, eps_r = 10.0),
        frequencies = 10.0 .^ range(log10(50.0), log10(1000.0); length = 101))
end

function retain(name, input, result)
    JLD2.jldsave(joinpath(OUTPUT, name * ".jld2"); problem = input, result)
    println(name, ": ", details(result).execution.source_run,
        " (reused=", details(result).execution.reused, ")")
    flush(stdout)
    @test size(result.Z) == size(result.Y) == (4, 4, 101)
    @test all(isfinite, result.Z.values) && all(isfinite, result.Y.values)
    @test result.f == input.frequencies
    settings = details(result).native_setting
    expected = Dict(string(component) => Dict(string(field) => control.readback
        for (field, control) in pairs(getproperty(settings, component)))
        for component in (:ground, :frequency, :configuration))
    @test details(result).native_readback == expected
    return result
end

println("Retained live checks: ", OUTPUT)
flush(stdout)
@testset "Live PSCAD normalization contracts" begin
    hot_problem = problem()
    labels = ("default", "default-tuple", "carson-pollaczek-lucca",
        "gary-wedepohl-ametani", "gary-saad-lucca")
    selections = [
        Formulation(:pscad),
        Formulation(:pscad; earth_impedance = (air = :default, earth = :default, mixed = :default)),
        Formulation(:pscad; earth_impedance = (
            air = :Carson1926, earth = :Pollaczek1926, mixed = :Lucca1994)),
        Formulation(:pscad; earth_impedance = (
            air = :Gary1976, earth = :WedepohlWilcox1973, mixed = :Ametani2009)),
        Formulation(:pscad; earth_impedance = (
            air = :Gary1976, earth = :Saad1996, mixed = :Lucca1994)),
    ]
    results = compute(hot_problem, selections; options = merge(OPTIONS,
        (on_result = (input, index, value) -> retain(labels[index], input, value),)))
    hot = first(results)
    cases = Set((row.kind, row.source, row.target)
        for row in details(hot).native_setting.interactions.earth_impedance)
    @test cases == Set(((:self, 1, 1), (:mutual, 1, 1), (:self, 2, 2),
        (:mutual, 2, 2), (:mutual, 1, 2), (:mutual, 2, 1)))
    for index in (2, 3)
        @test details(results[index]).execution.reused
        @test results[index].Z.values == hot.Z.values
        @test results[index].Z.values !== hot.Z.values
        @test results[index].Y.values == hot.Y.values
        @test details(results[index]).formulations.requested.earth_impedance ==
            LineCableModels.computation_details(selections[index]).requested.earth_impedance
    end

    lossy_selection = Formulation(:pscad; insulation_admittance = :Ametani2004)
    lossy = retain("lossy", hot_problem, compute(hot_problem, lossy_selection; options = OPTIONS))
    conductance = E.compare(hot, lossy, G)
    admittance = E.compare(hot, lossy, Y)
    @test all(ismissing, conductance.relative)
    @test all(==(:reference_below_tolerance), conductance.details.status)
    @test conductance.absolute[3, 3] > maximum(conductance.details.atol)
    @test all(!ismissing(admittance.relative[index, index]) for index in 1:4)
    @test all(ismissing, admittance.relative[1:2, 3:4])
    @test all(ismissing, admittance.relative[3:4, 1:2])

    screen_problem = problem(; semicon = true)
    screen_labels = ("semicon-lossless", "semicon-lossy")
    screens = compute(screen_problem,
        [Formulation(:pscad), Formulation(:pscad; semicon_admittance = :Ametani2004)];
        options = merge(OPTIONS, (on_result = (input, index, value) ->
            retain(screen_labels[index], input, value),)))
    screen_conductance = E.compare(screens[1], screens[2], G)
    @test screen_conductance.absolute[3, 3] > maximum(screen_conductance.details.atol)
    # Lossless constitutive inputs do not promise exactly zero native matrix G.
    # Preserve small native terms; numerical-zero treatment belongs to compare.
    for result in (hot, screens[1])
        document = EzXML.parsexml(details(result).exported_project)
        tangents = EzXML.findall("//User[@defn='master:Cable_Coax']/paramlist/param[@name='LT1']", document)
        @test length(tangents) == 4
        @test all(node -> iszero(parse(Float64, node["value"])), tangents)
    end

    cold_problem = problem(; temperature = 20)
    cold = retain("cold", cold_problem, compute(cold_problem, first(selections); options = OPTIONS))
    @test all(real(hot.Z.values[index, index, 1]) > real(cold.Z.values[index, index, 1]) for index in 1:4)
    # A constitutive temperature change and its explicitly supplied resistivity
    # produce identical native physical inputs. Force a separate solve to check
    # this equivalence through the native execution and import boundaries.
    equivalent_problem = problem(; temperature = 20, rho = 1.72e-8 * (1 + 0.004 * 40))
    equivalent = retain("equivalent-resistivity", equivalent_problem,
        compute(equivalent_problem, first(selections);
            options = merge(OPTIONS, (resume_run_directory = nothing,))))
    @test details(equivalent).exported_project == details(hot).exported_project
    @test equivalent.Z.values ≈ hot.Z.values rtol = 1e-10 atol = 0
    @test equivalent.Y.values ≈ hot.Y.values rtol = 1e-10 atol = 0

    recovered = retain("recovered", hot_problem,
        compute(hot_problem, first(selections); options = merge(OPTIONS,
            (resume_run_directory = details(hot).execution.source_run,))))
    @test details(recovered).execution.reused
    @test recovered.Z.values == hot.Z.values && recovered.Y.values == hot.Y.values

    # Different native formulations and loss laws are observations, not accuracy gates.
    comparisons = Dict(labels[index] => E.compare(hot, results[index]) for index in 4:5)
    JLD2.jldsave(joinpath(OUTPUT, "comparisons.jld2"); comparisons, conductance, admittance, screen_conductance)
    for (name, comparison) in sort!(collect(comparisons); by = first)
        println(name, " relative Z RMS per term: ", comparison.Z.relative)
    end
    println("Lossy G diagonal at first/last frequency: ",
        [(real(lossy.Y.values[index, index, 1]), real(lossy.Y.values[index, index, end])) for index in 1:4])
    println("Cold/hot R diagonal at first frequency: ",
        [(real(cold.Z.values[index, index, 1]), real(hot.Z.values[index, index, 1])) for index in 1:4])
end
