@testitem "ModalAnalysis / Newton eigenpair corrections and same-sample fallback" tags=[:unit] setup=[ModalTrackingExamples] begin
    using LinearAlgebra
    import LineCableModels.ModalAnalysis as M
    for R in (Float32, Float64)
        T=Complex{R}
        phase, expected=ModalTrackingExamples.phase_scan(R)
        selected=ModalAnalysisFormulation(:wedehpol1996)
        modal=@inferred compute(ModalAnalysisProblem(phase), selected)
        @test eltype(Ti(modal))===T
        @test gamma(modal)≈expected[[3, 2, 1], :]
        @test frequencies(modal)==phase.f
        @test Z(transform(PhaseDomain, modal))≈Z(phase)
        @test Y(transform(PhaseDomain, modal))≈Y(phase)
        @test H(PropagationParameters(modal))≈exp.(-25gamma(modal))
        for k in eachindex(phase.f)
            @test Y(phase)[:, :, k]*Z(phase)[:, :, k]*Ti(modal)[:, :, k]≈Ti(modal)[:, :, k]*Diagonal(gamma(modal)[:, k] .^
                                                                                                     2)
            paired=Z(phase)[:, :, k]*Ti(modal)[:, :, k]
            @test Tv(modal)[:, :, k]≈paired ./ transpose(norm.(eachcol(paired)))
            @test norm.(eachcol(Ti(modal)[:, :, k]))≈ones(R, 3)
        end
        diagnostics=details(modal).data.modal.diagnostics
        @test all(isnothing, diagnostics.iterations[:, 1])
        @test all(isfinite, diagnostics.eigen_residual)
        work=M.ModalAnalysisWorkspace(phase, selected.formula).buffers.newton
        # The kernel consumes an n+1 complex system. Use an independent diagonal problem.
        matrix=Matrix(Diagonal(T[1 + 0.2im, 3 + 0.7im, 5 + 0.1im]))
        vector=T[1, 0.01im, 0.02]
        tolerance=R===Float32 ? R(1e-5) : R(1e-10)
        value, converged, iterations=M.newton_eigenpair!(vector, T(0.9+0.1im), matrix,
            (convergence = tolerance, max_iterations = 60), work)
        @test converged && iterations>0
        @test value≈matrix[1, 1]
        @test norm(matrix*vector-value*vector)<10tolerance
        @test norm(vector)≈one(R)
        _, converged, _=M.newton_eigenpair!(T[1, 0, 0], zero(T), zeros(T, 3, 3),
            (convergence = tolerance, max_iterations = 60), work)
        @test !converged
    end
    # Both Newton trajectories reach the same eigenpair. Direct fallback restores rank.
    z=cat(Matrix{ComplexF64}(I, 2, 2), ComplexF64[2 -1; -1 2]; dims = 3)
    y=cat(Matrix{ComplexF64}(I, 2, 2), Matrix{ComplexF64}(I, 2, 2); dims = 3)
    collapsed=compute(ModalAnalysisProblem(LineParameters(z, y, [1.0, 2.0])),
        ModalAnalysisFormulation(:wedehpol1996))
    diagnostics=details(collapsed).data.modal.diagnostics
    @test diagnostics.fallback_frequencies==[2]
    @test diagnostics.missed_frequencies==[2]
    @test all(diagnostics.converged[:, 2])
    @test sort(real.(gamma(collapsed)[:, 2]) .^ 2)≈[1, 3]
    @test rank(Ti(collapsed)[:, :, 2])==2
    phase, _=ModalTrackingExamples.phase_scan(; angle_step = 0.4)
    limited=compute(ModalAnalysisProblem(phase),
        ModalAnalysisFormulation(:wedehpol1996;
            options = (iteration = (max_iterations = 1,),)))
    @test !isempty(details(limited).data.modal.diagnostics.fallback_frequencies)
    for k in eachindex(phase.f)
        @test Y(phase)[:, :, k]*Z(phase)[:, :, k]*Ti(limited)[:, :, k]≈Ti(limited)[:, :, k]*Diagonal(gamma(limited)[
            :, k] .^ 2)
    end
    # A zero product is a valid normalized eigensystem. It must not divide by zero.
    zero_product=compute(
        ModalAnalysisProblem(LineParameters(ones(ComplexF64, 1, 1, 2),
            zeros(ComplexF64, 1, 1, 2), [1.0, 2.0])),
        ModalAnalysisFormulation(:wedehpol1996))
    @test iszero(gamma(zero_product))
    @test H(PropagationParameters(zero_product; line_length = 10.0))==ones(ComplexF64, 1, 2)
    for invalid in ((convergence = 0,), (convergence = NaN,), (max_iterations = 0,),
        (max_iterations = true,), (max_iterations = 1.5,), (unused = 1,))
        @test_throws ArgumentError ModalAnalysisFormulation(:wedehpol1996; options = (iteration = invalid,))
    end
    @test formula_id(ModalAnalysisFormulation().formula)===:chrysochos2014
    @test :wedehpol1996 in M.formulas()
end

@testitem "ModalAnalysis / shared paired rotation and sign continuity" tags=[:unit] begin
    using LinearAlgebra
    import LineCableModels.ModalAnalysis as M
    import LineCableModels.Engine as E
    import LineCableModels.Grammar as G
    # A 60-degree scalar catches the supplied routine's second application of the angle.
    for R in (Float32, Float64)
        T=Complex{R}
        current=reshape(T[cis(R(pi/3))], 1, 1, 1)
        voltage=2im .* current
        M.orient_modes!(voltage, current, true)
        @test real(only(current))≈one(R)
        @test abs(imag(only(current)))<8eps(R)
        @test only(voltage)≈2im
        @test eltype(current)===T
        for exponent in (-1, 1)
            scale=R(10)^(exponent*(R===Float32 ? 25 : 200))
            current=reshape(T[scale * cis(R(pi / 3)), -scale * cis(R(pi / 3))], 1, 1, 2)
            voltage=copy(current)
            M.orient_modes!(voltage, current, true)
            @test vec(current) ./ scale ≈ ones(T, 2)
            @test all(isfinite, current)
        end
        original=reshape(T[1 + 2im, 3 - 0.4im, 0.2 + im, 2 - 3im], 2, 2, 1)
        current=copy(original)
        voltage=3 .* original
        M.orient_modes!(voltage, current, true)
        for mode in 1:2
            before=original[:, mode, 1]
            after=current[:, mode, 1]
            minimum_imaginary=(sum(abs2, before)-abs(sum(before .^ 2)))/2
            @test sum(abs2, imag.(after))≈minimum_imaginary rtol=32eps(R)
            @test norm(after)≈norm(before)
        end
        @test voltage≈3current
    end
    struct FixedOrientation <: AbstractFormulation
        voltage::Array{ComplexF64, 3}
        current::Array{ComplexF64, 3}
        roots::Matrix{ComplexF64}
    end
    LineCableModels.formula_id(::FixedOrientation) = :fixed_orientation
    Base.NamedTuple(::FixedOrientation) = (identifier = :fixed_orientation,)
    G.formulation_options(::FixedOrientation) = FormulationOptions()
    E.initialize_buffers(::FixedOrientation, ::Type, input, invariants, buffers) = buffers
    function M.decompose!(selected::FixedOrientation, workspace, parameters, options)
        copyto!(workspace.Tv, selected.voltage)
        copyto!(workspace.Ti, selected.current)
        copyto!(workspace.roots, selected.roots)
        return workspace
    end
    f=[10.0, 20.0, 30.0]
    C=ComplexF64[1+0.7im 0.3-0.2im; 0.2+0.3im 1-0.4im]
    z=ComplexF64[2 + 3im, 4 + im]
    y=ComplexF64[0.1 + 0.2im, 0.2 + 0.3im]
    impedance=C*Diagonal(z)/C
    admittance=C*Diagonal(y)/C
    current=cat(C, -C, cis(0.2) .* C; dims = 3)
    voltage=copy(current)
    roots=repeat(sqrt.(z .* y), 1, 3)
    phase=LineParameters(
        cat(fill(impedance, 3)...; dims = 3), cat(fill(admittance, 3)...; dims = 3), f;
        details = ComputationDetails((
            gridpoint = G.gridpoint_id(), selections = (Z = :example, Y = :example))))
    selected=ModalAnalysisFormulation(FixedOrientation(voltage, current, roots))
    rotated=@inferred compute(ModalAnalysisProblem(phase), selected)
    raw=compute(ModalAnalysisProblem(phase), selected; options = (rotate = false,))
    @test details(rotated).data.modal.rotate===true
    @test details(raw).data.modal.rotate===false
    @test Ti(raw)[:, :, 1]==C
    @test Ti(raw)[:, :, 2]==C # the 180-degree flip is removed even with rotation disabled
    @test Ti(rotated)[:, :, 1]≈Ti(rotated)[:, :, 2]≈Ti(rotated)[:, :, 3]
    @test gamma(rotated)==gamma(raw)==roots
    @test Zc(rotated)≈Zc(raw)
    @test Yc(rotated)≈Yc(raw)
    for result in (rotated, raw)
        @test Z(transform(PhaseDomain, result))≈Z(phase)
        @test Y(transform(PhaseDomain, result))≈Y(phase)
        for k in 2:3, mode in 1:2

            @test real(dot(Ti(result)[:, mode, k - 1], Ti(result)[:, mode, k]))>=0
        end
    end
    @test selected.formula.current==current
    first_segment=PropagationParameters(rotated; line_length = 1.0)
    second_segment=PropagationParameters(raw; line_length = 1.0)
    for field in (:current, :voltage)
        @test H(first_segment, PhaseDomain; field)≈H(second_segment, PhaseDomain; field)
    end
    @test Zc(rotated, PhaseDomain)≈Zc(raw, PhaseDomain)
    @test Yc(rotated, PhaseDomain)≈Yc(raw, PhaseDomain)
    observed=observables.([first_segment, second_segment], Ref((Ti,)))
    @test length(G.observation_groups(observed; request = (Ti, real)))==2
    mktempdir() do folder
        for extension in ("json", "jls")
            path=joinpath(folder, "rotation."*extension)
            LineCableModels.save(observed[1], path)
            restored=LineCableModels.import_data(Val(:observed), path)
            @test observe(restored, (Ti, real))==observe(observed[1], (Ti, real))
            @test G.observation_product(restored, (Ti, real)).assumptions.rotate===true
        end
    end
    for invalid in (0, 1, :yes, nothing)
        @test_throws ArgumentError compute(ModalAnalysisProblem(phase), selected; options = (rotate = invalid,))
    end
end

@testitem "ModalAnalysis / Newton collection, convenience and uncertain consumers" tags=[:unit] setup=[
    ModalTrackingExamples, TestFixtures] begin
    using Measurements, LinearAlgebra
    phase, _=ModalTrackingExamples.phase_scan()
    selected=ModalAnalysisFormulation(:wedehpol1996)
    source=ParametricResult(Combinatorial(selected), [phase, phase])
    results=compute(Gridspace{ModalAnalysisProblem}(source), selected; options = (rotate = false,))
    @test length(results)==2 && isconcretetype(eltype(results))
    @test all(result->details(result).data.modal.rotate===false, results)
    @test gamma(results[1])==gamma(results[2])
    total=LineParameters(PhaseDomain, SeriesImpedance(Z(phase) .* 25; basis = :total),
        ShuntAdmittance(Y(phase) .* 25; basis = :total), phase.f, phase.details)
    @test H(PropagationParameters(compute(ModalAnalysisProblem(total), selected)))≈H(PropagationParameters(results[1]))
    retained=results[1][[4, 2]]
    @test gamma(retained)==gamma(results[1])[:, [4, 2]]
    @test details(retained).data.modal.rotate===false
    observed=observables(PropagationParameters(retained), (gamma, Zc, H, Tv, Ti))
    @test size(report(observed)[alpha])==(2, 4)
    x=measurement(2.0, 0.1)
    uncertain=LineParameters(reshape([complex(x, 0)], 1, 1, 1), reshape([3+0im], 1, 1, 1), [50.0])
    modal=@inferred compute(ModalAnalysisProblem(uncertain), selected)
    @test Measurements.derivative(real(only(gamma(modal))), x)≈3/(2sqrt(6))
    problem=TestFixtures.three_bare_wires_problem(frequencies = [50.0, 80.0])
    phase=compute(problem, Formulation())
    direct=compute(ModalAnalysisProblem(phase), selected; options = (rotate = false,))
    upstream=Ref(0)
    downstream=Ref(0)
    combined=compute(problem, Formulation(); modal = :wedehpol1996,
        options = (on_result = (args...)->(upstream[]+=1),),
        modal_options = (rotate = false, on_result = (args...)->(downstream[]+=1)))
    @test upstream[]==downstream[]==1
    @test gamma(combined)==gamma(direct)
    @test Ti(combined)==Ti(direct)
    @test details(combined).data.modal.rotate===false
end
