@testset "Engine matrices inside the cables" begin
    Engine = LineCableModels.Engine
    reductions = LineCableModels.FormulationOptions(CoaxialFormulation)
    internal = Engine.InternalImpedance.Formula(Val(:schelkunoff1934))
    term = Engine.InsulationImpedance.Formula(Val(:ametani1980))
    shunt = Engine.Shunt.Formula(Val(:series))
    law = Engine.Resistivity.Formula(Val(:linear))
    insulating = Engine.InsulationAdmittance.Formula(Val(:lossless))
    semiconducting = Engine.SemiconAdmittance.Formula(Val(:lossless))

    @testset "$case against its stage targets" for (case, problem) in REFERENCES
        input = Engine.CoaxialInput(problem, reductions)
        n = length(input.cable)
        c = length(input.radius)
        targets = Dict(matrix => zeros(ComplexF64, rows, rows, length(problem.frequencies))
        for (matrix, rows) in
            (("Zin", n), ("Pin", n), ("Zg", c), ("Pg", c), ("Z", n), ("Y", n)))
        for line in eachline(joinpath(@__DIR__, "..", "references", case, "stages.csv"))
            fields = split(line, ',')
            startswith(line, "analytical") && haskey(targets, fields[3]) || continue
            k = searchsortedfirst(problem.frequencies, parse(Float64, fields[2]))
            targets[fields[3]][parse(Int, fields[4]), parse(
                Int, fields[5]), k] = complex(
                parse(Float64, fields[6]), parse(Float64, fields[7]))
        end
        for line in eachline(joinpath(@__DIR__, "..", "references", case, "matrices.csv"))
            fields = split(line, ',')
            startswith(line, "analytical") || continue
            k = searchsortedfirst(problem.frequencies, parse(Float64, fields[2]))
            targets[fields[3]][parse(Int, fields[4]), parse(
                Int, fields[5]), k] = complex(
                parse(Float64, fields[6]), parse(Float64, fields[7]))
        end
        rho = [law(tube, problem.temperature)
               for tube in input.material]
        semicons = findall(layer -> layer isa Material{Float64, :semicon},
            input.layer_material)
        Zin = Array{ComplexF64}(undef, n, n, length(problem.frequencies))
        Pin = similar(Zin)
        for (k, f) in enumerate(problem.frequencies)
            s = im * (2π * f)
            kappa = ComplexF64[(i in semicons ? semiconducting : insulating)(layer, f)
                               for (i, layer) in enumerate(input.layer_material)]
            Engine.impedance!(internal, term, view(Zin,:,:,k), input, rho, s)
            Engine.potential!(shunt, view(Pin,:,:,k), input, kappa, s)
        end
        for (computed, target) in ((Zin, targets["Zin"]), (Pin, targets["Pin"]))
            @test all(maximum(abs, computed[:, :, k] - target[:, :, k]) <=
                      1.0e-12 * maximum(abs(target[i, i, k]) for i in 1:n)
            for k in eachindex(problem.frequencies))
        end
        result = compute(problem, CoaxialFormulation(); trace = true)
        for (matrix, computed, tolerance) in (("Z", result.Z, 1e-9), ("Y", result.Y, 1e-9),
            ("Zin", result.trace.Zin, 1e-12), ("Pin", result.trace.Pin, 1e-12),
            ("Zg", result.trace.Zg, 1e-9), ("Pg", result.trace.Pg, 1e-9))
            target = targets[matrix]
            errors = [maximum(abs, computed[:, :, k] - target[:, :, k])
                      for k in eachindex(problem.frequencies)]
            scales = [maximum(abs(target[i, i, k]) for i in axes(target, 1))
                      for k in eachindex(problem.frequencies)]
            @test all(errors .<= tolerance .* scales)
        end
    end
end
