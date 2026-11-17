@testset "Engine matrices inside the cables" begin
    Engine = LineCableModels.Engine
    reductions = LineCableModels.FormulationOptions(CoaxialFormulation)
    internal = Engine.InternalImpedance.Formula(Val(:schelkunoff1934))
    term = Engine.InsulationImpedance.Formula(Val(:ametani1980))
    shunt = Engine.Shunt.Formula(Val(:series))
    law = Engine.Resistivity.Formula(Val(:linear))
    insulating = Engine.InsulationAdmittance.Formula(Val(:lossless))
    semiconducting = Engine.SemiconAdmittance.Formula(Val(:lossless))
    # The order of the formulation of a case's analytical rows, where it is not the default.
    orders = Dict("three_cores_steel_pipe" => 20)

    @testset "$case against its rows of record" for (case, problem) in REFERENCES
        result = compute(
            problem, CoaxialFormulation(order = get(orders, case, 0)); trace = true)
        computed = Dict("Z" => result.Z, "Y" => result.Y, "Zin" => result.trace.Zin,
            "Pin" => result.trace.Pin, "Zg" => result.trace.Zg, "Pg" => result.trace.Pg)
        # The file `matrices.csv` of a case gives the rows of Z and Y. Its `stages.csv`,
        # where the folder has one, gives the rows of Zin, Pin, Zg and Pg.
        stages = joinpath(@__DIR__, "..", "references", case, "stages.csv")
        targets = Dict(matrix => zeros(ComplexF64, size(computed[matrix]))
        for matrix in (isfile(stages) ? ("Zin", "Pin", "Zg", "Pg", "Z", "Y") : ("Z", "Y")))
        for file in (stages, joinpath(dirname(stages), "matrices.csv"))
            isfile(file) || continue
            for line in eachline(file)
                fields = split(line, ',')
                startswith(line, "analytical") && haskey(targets, fields[3]) || continue
                k = searchsortedfirst(problem.frequencies, parse(Float64, fields[2]))
                targets[fields[3]][parse(Int, fields[4]), parse(
                    Int, fields[5]), k] = complex(
                    parse(Float64, fields[6]), parse(Float64, fields[7]))
            end
        end
        # Each matrix lies within 1e-12 of its rows of record, relative to the largest
        # diagonal entry of the record at the frequency. The rectangular core's Zin rows
        # use gmd at τ = 1e-6, measured on Julia 1.12.7, x86_64 Linux.
        @testset "$matrix" for (matrix, target) in targets
            @test all(maximum(abs, computed[matrix][:, :, k] - target[:, :, k]) <=
                      1.0e-12 * maximum(abs(target[i, i, k]) for i in axes(target, 1))
            for k in eachindex(problem.frequencies))
        end
        # The recipes of the formulation, called with the default formulas, give the rows
        # of Zin and Pin within the same bound.
        if isfile(stages)
            input = Engine.CoaxialInput(problem, reductions)
            buffers = Engine.initialize_buffers(CoaxialFormulation(), problem, input)
            n = length(input.cable)
            rho = [law(tube, problem.temperature)
                   for tube in input.material]
            semicons = findall(dielectric -> dielectric isa Material{Float64, :semicon},
                input.dielectric_material)
            Zin = Array{ComplexF64}(undef, n, n, length(problem.frequencies))
            Pin = similar(Zin)
            for (k, f) in enumerate(problem.frequencies)
                s = im * (2π * f)
                kappa = ComplexF64[(i in semicons ? semiconducting : insulating)(dielectric, f)
                                   for (i, dielectric) in
                                       enumerate(input.dielectric_material)]
                Engine.impedance!(internal, term, view(Zin,:,:,k), buffers, input, rho, s)
                Engine.potential!(shunt, view(Pin,:,:,k), buffers, input, kappa, s)
            end
            for (matrix, target) in ((Zin, targets["Zin"]), (Pin, targets["Pin"]))
                @test all(maximum(abs, matrix[:, :, k] - target[:, :, k]) <=
                          1.0e-12 * maximum(abs(target[i, i, k]) for i in 1:n)
                for k in eachindex(problem.frequencies))
            end
        end
    end
end
