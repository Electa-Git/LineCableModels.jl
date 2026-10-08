@testmodule ModalTrackingExamples begin
    using LineCableModels, LinearAlgebra

    function phase_scan(::Type{R} = Float64; clustered = false, angle_step = 0.02) where {R}
        f = R[50, 80, 130, 200]
        Z = zeros(Complex{R}, 3, 3, length(f))
        Y = similar(Z)
        roots = Matrix{Complex{R}}(undef, 3, length(f))
        for k in eachindex(f)
            angle = R(angle_step) * (k - 1)
            rotation = R[cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)]
            z = Complex{R}[2 + im * f[k] / 20,
                (clustered ? 2 : 4) + im * f[k] / (clustered ? 20 : 10),
                8 + im * f[k] / 5]
            y = complex(R(1e-7), R(1e-8) * f[k])
            Z[:, :, k] = rotation * Diagonal(z) * transpose(rotation)
            Y[:, :, k] = y * Matrix{R}(I, 3, 3)
            roots[:, k] = sqrt.(z .* y)
        end
        phase = LineParameters(Z, Y, f;
            details = ComputationDetails((inputs = (system = (line_length = R(25),),),)))
        return phase, roots
    end
end

@testitem "ModalAnalysis / Vieira complex residual, Jacobian and minimum-norm steps" tags=[:unit, :modal] begin
    using LinearAlgebra
    import LineCableModels.ModalAnalysis as M
    formula = M.Formula(:vieira2026)
    for R in (Float32, Float64)
        T = Complex{R}
        S = T[1+im 2-im; 0.5im 3+2im]
        x = T[0.8 + 0.2im, 0.3 - 0.1im, 1.2 + 0.4im]
        residual = zeros(T, 3)
        J = zeros(T, 3, 3)
        @test M.eigenpair_residual!(residual, x, S) === residual
        @test residual ≈ vcat((S-x[3]*I)*x[1:2], sum(x[1:2] .^ 2)-1)
        M.eigenpair_jacobian!(J, x, S)
        expected = T[S[1, 1]-x[3] S[1, 2] -x[1];
                     S[2, 1] S[2, 2]-x[3] -x[2]; 2x[1] 2x[2] 0]
        @test J == expected
        if R === Float64
            direction = T[0.1 + 0.3im, 0.2 - 0.2im, -0.4 + 0.1im]
            plus, minus = similar(residual), similar(residual)
            M.eigenpair_residual!(plus, x+1e-5direction, S)
            M.eigenpair_residual!(minus, x-1e-5direction, S)
            @test (plus-minus)/2e-5 ≈ J*direction rtol=1e-9
        end
        A = T[1 im; 2 2im; 0 0]
        b = T[3 + im, 6 + 2im, 4]
        solution, projection = zeros(T, 2), zeros(T, 2)
        M.minimum_norm!(solution, copy(A), b, projection)
        @test solution ≈ T[(3 + im) / 2, (1 - 3im) / 2]
        coefficients = zeros(T, 2, 2)
        M.minimum_norm!(coefficients, copy(A), hcat(b, 2b), similar(coefficients))
        @test coefficients ≈ hcat(solution, 2solution)
        M.minimum_norm!(solution, zeros(T, 3, 2), b, projection)
        @test iszero(solution)

        least_squares = LineCableModels.Commons.initialize_buffers(
            formula, T, (;), (; n = 2, nf = 1), (;)).least_squares
        vector = T[1, 0.1im]
        matrix = Matrix(Diagonal(T[1 + 0.2im, 3 + 0.7im]))
        tolerance = R === Float64 ? 1e-11 : 1e-5
        value, converged, iterations = M.levenberg_marquardt_step!(formula,
            vector, T(0.9+0.1im), matrix,
            (convergence = tolerance, max_iterations = 100), least_squares)
        @test converged
        @test iterations > 0
        @test value ≈ matrix[1, 1]
        @test sum(vector .^ 2) ≈ one(T)
        @test norm(matrix*vector-value*vector) < tolerance
        # An exact eigenpair does not need iterations. A stationary non-root stalls.
        value, converged, iterations = M.levenberg_marquardt_step!(formula,
            T[1, 0], matrix[1, 1], matrix, (convergence = tolerance, max_iterations = 100), least_squares)
        @test converged && iterations == 0
        vector = zeros(T, 2)
        _, converged, iterations = M.levenberg_marquardt_step!(formula,
            vector, zero(T), matrix, (convergence = tolerance, max_iterations = 2), least_squares)
        @test !converged && iterations == 1
        @test iszero(vector)
    end
    cost = [1.0 2; 2 100]
    assignment = zeros(Int, 2)
    @test M.greedy_assignment!(assignment, cost) == [1, 2]
end

@testitem "ModalAnalysis / Vieira tracked modes, ordering and clustered subspaces" tags=[:unit, :modal, :slow] setup=[ModalTrackingExamples] begin
    using LinearAlgebra
    for R in (Float32, Float64), clustered in (false, true)

        phase, expected_roots=ModalTrackingExamples.phase_scan(R; clustered)
        before_Z, before_Y=copy(Z(phase)), copy(Y(phase))
        modal=@inferred compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:vieira2026))
        @test eltype(Ti(modal)) === Complex{R}
        @test eltype(gamma(modal)) === Complex{R}
        @test gamma(modal) ≈ expected_roots[[3, 2, 1], :]
        @test all(>=(zero(R)), real.(gamma(modal)))
        for k in eachindex(phase.f)
            current, voltage=Ti(modal)[:, :, k], Tv(modal)[:, :, k]
            @test norm.(eachcol(current)) ≈ ones(R, 3)
            @test norm.(eachcol(voltage)) ≈ ones(R, 3)
            @test Y(phase)[:, :, k]*Z(phase)[:, :, k]*current ≈
                  current*Diagonal(gamma(modal)[:, k] .^ 2)
            paired=Z(phase)[:, :, k]*current
            @test voltage ≈ paired ./ norm.(eachcol(paired))'
            if clustered
                angle=R(0.02)*(k-1)
                target=R[cos(angle) 0; 0 1; -sin(angle) 0]
                repeated=current[:, 2:3]
                @test repeated*pinv(repeated) ≈ target*transpose(target)
            end
        end
        rebuilt=transform(PhaseDomain, modal)
        @test Z(rebuilt) ≈ before_Z
        @test Y(rebuilt) ≈ before_Y
        @test Z(phase) == before_Z && Y(phase) == before_Y
        diagnostics=details(modal).data.modal.diagnostics
        @test all(isfinite, diagnostics.eigen_residual)
        @test all(isnothing, diagnostics.iterations[:, 1])
        @test all(x -> x isa Int && x >= 0, diagnostics.iterations[:, 2:end])
        unordered=compute(ModalAnalysisProblem(phase),
            ModalAnalysisFormulation(:vieira2026;
                options = (tracking = (order_by_velocity = false,),)))
        order=sortperm(imag.(gamma(unordered)[:, end]); rev = true)
        @test gamma(modal) ≈ gamma(unordered)[order, :]
        @test Ti(modal) ≈ Ti(unordered)[:, order, :]
        @test isequal(diagnostics.iterations, details(unordered).data.modal.diagnostics.iterations[order, :])
        @test isequal(diagnostics.converged, details(unordered).data.modal.diagnostics.converged[order, :])
        @test diagnostics.eigen_residual ≈
              details(unordered).data.modal.diagnostics.eigen_residual[order, :]
    end
end

@testitem "ModalAnalysis / Vieira same-sample fallback and numerical controls" tags=[:unit, :modal] setup=[ModalTrackingExamples] begin
    using LinearAlgebra
    import LineCableModels.ModalAnalysis as M
    phase, roots=ModalTrackingExamples.phase_scan(; angle_step = 0.4)
    selected=ModalAnalysisFormulation(:vieira2026;
        options = (
            iteration = (max_iterations = 1,), tracking = (predictor_tolerance = 1e-8,)))
    modal=compute(ModalAnalysisProblem(phase), selected)
    diagnostics=details(modal).data.modal.diagnostics
    @test !isempty(diagnostics.fallback_frequencies)
    @test all(in(diagnostics.missed_frequencies), diagnostics.fallback_frequencies)
    @test modal.f == phase.f
    for k in eachindex(phase.f)
        @test Y(phase)[:, :, k]*Z(phase)[:, :, k]*Ti(modal)[:, :, k] ≈
              Ti(modal)[:, :, k]*Diagonal(gamma(modal)[:, k] .^ 2)
    end
    strict_match=compute(ModalAnalysisProblem(phase),
        ModalAnalysisFormulation(:vieira2026;
            options = (iteration = (max_iterations = 1,),
                tracking = (eigenvalue_tolerance = 1e-15,))))
    @test !isempty(details(strict_match).data.modal.diagnostics.fallback_frequencies)
    @test formula_id(ModalAnalysisFormulation().formula) === :chrysochos2014
    @test :vieira2026 in M.formulas(M.Formula)
    for controls in ((convergence = 0,), (convergence = NaN,), (max_iterations = 0,),
        (max_iterations = true,), (max_iterations = 1.5,), (unused = 1,))
        @test_throws ArgumentError ModalAnalysisFormulation(:vieira2026; options = (iteration = controls,))
    end
    for controls in ((predictor_tolerance = 0,), (eigenvalue_tolerance = Inf,),
        (cluster_tolerance = -1,), (order_by_velocity = 1,), (max_depth = 12,))
        @test_throws ArgumentError ModalAnalysisFormulation(:vieira2026; options = (tracking = controls,))
    end
    @test_throws ArgumentError ModalAnalysisFormulation(:vieira2026; parameters = (s = [1im],))
    branch=LineParameters(reshape(ComplexF64[-im], 1, 1, 1),
        reshape(ComplexF64[-im], 1, 1, 1), [50.0])
    @test only(gamma(compute(ModalAnalysisProblem(branch), ModalAnalysisFormulation(:vieira2026)))) ≈
          im
end

@testitem "ModalAnalysis / Vieira public selection, propagation, observations and uncertainty" tags=[:unit, :measurements, :slow] setup=[ModalTrackingExamples] begin
    using Measurements, LinearAlgebra
    phase, roots=ModalTrackingExamples.phase_scan()
    selected=ModalAnalysisFormulation(:vieira2026)
    modal=compute(ModalAnalysisProblem(phase), selected)
    segment=PropagationParameters(modal)
    @test H(segment) ≈ exp.(-gamma(modal) .* 25)
    total=LineParameters(PhaseDomain, SeriesImpedance(Z(phase) .* 25; basis = :total),
        ShuntAdmittance(Y(phase) .* 25; basis = :total), phase.f, phase.details)
    total_modal=compute(ModalAnalysisProblem(total), selected)
    @test gamma(total_modal) ≈ gamma(modal) .* 25
    @test H(PropagationParameters(total_modal)) ≈ H(segment)
    @test Zc(total_modal) ≈ Zc(modal)
    for selection in (2, 2:3, [4, 2], :)
        indices=selection isa Int ? (selection:selection) : selection
        retained=modal[selection]
        @test gamma(retained) == gamma(modal)[:, indices]
        @test Ti(retained) == Ti(modal)[:, :, indices]
        @test details(retained).data.modal.diagnostics.eigen_residual ==
              details(modal).data.modal.diagnostics.eigen_residual[:, indices]
        @test H(segment[selection]) == H(segment)[:, indices]
    end
    observed=observables(segment, (gamma, Zc, H, Tv, Ti))
    @test observe(observed, (gamma, real)) ≈ 1000 .* observe(segment, gamma, real)
    tables=report(observed)
    @test size(tables[(gamma, real)]) == (length(phase.f), 4)
    @test tables[alpha] === tables[(gamma, real)]
    mktempdir() do folder
        for extension in ("json", "jls")
            path=joinpath(folder, "vieira."*extension)
            LineCableModels.save(observed, path)
            restored=LineCableModels.import_data(Val(:observed), path)
            @test observe(restored, (H, real)) == observe(observed, (H, real))
        end
    end
    x=measurement(2.0, 0.1)
    uncertain=LineParameters(reshape([complex(x, zero(x))], 1, 1, 1),
        reshape([3.0+0im], 1, 1, 1), [50.0])
    propagated=@inferred compute(ModalAnalysisProblem(uncertain), selected)
    root=real(only(gamma(propagated)))
    @test nominal(root) ≈ sqrt(6)
    @test Measurements.derivative(root, x) ≈ 3/(2sqrt(6))
    factor=real(only(H(PropagationParameters(propagated; line_length = 10.0))))
    @test Measurements.derivative(factor, x) ≈ -10exp(-10sqrt(6))*3/(2sqrt(6))
    source=ParametricResult(Combinatorial(selected), [phase, phase])
    results=compute(Gridspace{ModalAnalysisProblem}(source), selected)
    @test length(results) == 2
    @test isconcretetype(eltype(results))
    @test gamma.(results) == [gamma(modal), gamma(modal)]
    variants=ModalAnalysisFormulation(Grid((:chrysochos2014, :vieira2026)))
    combinations=compute(Gridspace{ModalAnalysisProblem}(source), variants)
    @test length(combinations) == 4
    @test [details(result).data.modal.identifier for result in combinations] ==
          [:chrysochos2014, :chrysochos2014, :vieira2026, :vieira2026]
end

@testitem "ModalAnalysis / Vieira convenience computation preserves completed upstream results" tags=[:unit, :parametric] setup=[TestFixtures] begin
    problem=TestFixtures.three_bare_wires_problem(frequencies = [50.0, 80.0])
    formulation=Formulation()
    phase=compute(problem, formulation)
    direct=compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:vieira2026))
    upstream=Ref(0)
    downstream=Ref(0)
    combined=compute(problem, formulation; modal = :vieira2026,
        options = (trace = true, on_result = (args...)->(upstream[]+=1)),
        modal_options = (on_result = (args...)->(downstream[]+=1),))
    @test upstream[] == downstream[] == 1
    @test gamma(combined) == gamma(direct)
    @test Ti(combined) == Ti(direct)
    @test Z(combined) == Z(direct)
    @test haskey(details(combined).data, :trace)
    @test !haskey(details(direct).data, :trace)
end
