@testitem "Transforms / repeated eigenvalues preserve the modal subspace" tags=[:unit] begin
    using LinearAlgebra

    frequencies = [50.0, 55.0, 60.0, 65.0, 70.0]
    impedance = zeros(ComplexF64, 3, 3, length(frequencies))
    admittance = similar(impedance)
    principal_axes = Matrix{Float64}[]
    for (index, angle) in enumerate(range(0.0, 0.3; length = length(frequencies)))
        rotation = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)]
        push!(principal_axes, rotation)
        impedance[:, :, index] = (1 + 2im) .* (
            rotation * Diagonal([1.0, 1.0, 5.0]) * transpose(rotation))
        admittance[:, :, index] = (1e-8 + 4e-7im) .* Matrix(I, 3, 3)
    end
    phase = LineParameters(impedance, admittance, frequencies)
    modal = @inferred compute(ModalTransformationProblem(phase),
        ModalTransformationFormulation(:default))
    maps = operators(modal)
    for index in eachindex(frequencies)
        vectors = transpose(maps.voltage[:, :, index])
        product = admittance[:, :, index] * impedance[:, :, index]
        eigenvalues = diag(maps.current[:, :, index] * product * vectors)
        @test norm(product * vectors - vectors * Diagonal(eigenvalues)) <=
              1e-8 * norm(product) * norm(vectors)
        @test maps.current[:, :, index] * vectors ≈ Matrix(I, 3, 3) rtol=1e-12
        # Only the two-dimensional repeated eigenspace is unique, not its basis.
        target = (1e-8 + 4e-7im) * (1 + 2im)
        repeated = findall(value -> abs(value - target) <= 1e-8 * abs(target), eigenvalues)
        @test length(repeated) == 2
        actual_basis = vectors[:, repeated]
        actual_projector = actual_basis * pinv(actual_basis)
        expected_basis = principal_axes[index][:, 1:2]
        @test actual_projector ≈ expected_basis * transpose(expected_basis) rtol=1e-8
    end
    rebuilt = @inferred compute(ModalTransformationProblem(modal))
    @test rebuilt.Z.values ≈ impedance rtol=1e-12
    @test rebuilt.Y.values ≈ admittance rtol=1e-12
    @test phase.Z.values == impedance
    @test phase.Y.values == admittance

    # Frequency selections retain the original modal basis at each selected sample.
    for selector in (2, 2:4, [5, 2], :)
        indices = selector isa Integer ? (selector:selector) : selector
        selected = @inferred modal[selector]
        selected_maps = operators(selected)
        @test selected.f == frequencies[indices]
        @test selected.Z.values == modal.Z.values[:, :, indices]
        @test selected.Y.values == modal.Y.values[:, :, indices]
        @test selected_maps.voltage == maps.voltage[:, :, indices]
        @test selected_maps.current == maps.current[:, :, indices]
        selected_phase = @inferred compute(ModalTransformationProblem(selected))
        @test selected_phase.Z.values ≈ impedance[:, :, indices] rtol=1e-12
        @test selected_phase.Y.values ≈ admittance[:, :, indices] rtol=1e-12
    end

    selected = modal[2:4]
    voltage_before = copy(maps.voltage)
    current_before = copy(maps.current)
    operators(selected).voltage[1, 1, 1] += 1
    operators(selected).current[1, 1, 1] += 1
    @test maps.voltage == voltage_before
    @test maps.current == current_before
end

@testitem "Transforms / limited iteration retains an explicit matched eigensolution" tags=[:unit] begin
    using LinearAlgebra

    frequencies = [50.0, 100.0, 200.0, 400.0]
    impedance = zeros(ComplexF64, 2, 2, length(frequencies))
    admittance = similar(impedance)
    for (index, angle) in enumerate((0.0, 0.4, 0.8, 1.2))
        rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)]
        impedance[:, :, index] = rotation *
                                 Diagonal(ComplexF64[
                                 1 + index * im, 3 + 2index * im]) * transpose(rotation)
        admittance[:, :, index] = rotation *
                                  Diagonal(ComplexF64[
                                  1e-8 + index * 1e-7im, 2e-8 + index * 3e-7im]) *
                                  transpose(rotation)
    end
    phase = LineParameters(impedance, admittance, frequencies)
    for identifier in (:default,)
        @testset "$identifier" begin
            formulation = ModalTransformationFormulation(
                identifier; options = (iteration = (max_iterations = 1,),))
            modal = @test_logs (:warn, r"retained matched eigensolutions") match_mode=:any compute(
                ModalTransformationProblem(phase), formulation)
            @test !isempty(details(modal).modal.fallback_frequencies)
            @test details(modal).modal.options.iteration.fallback === :matched
            @test_throws ErrorException compute(ModalTransformationProblem(phase),
                ModalTransformationFormulation(identifier;
                    options =
                    (iteration = (max_iterations = 1, fallback = :error),)))
            maps = operators(modal)
            for index in eachindex(frequencies)
                vectors = transpose(maps.voltage[:, :, index])
                product = admittance[:, :, index] * impedance[:, :, index]
                eigenvalues = diag(maps.current[:, :, index] * product * vectors)
                @test norm(product * vectors - vectors * Diagonal(eigenvalues)) <=
                      1e-10 * norm(product) * norm(vectors)
                @test maps.current[:, :, index] * vectors ≈ Matrix(I, 2, 2) rtol=1e-12
                @test all(isfinite, modal.Z.values[:, :, index])
                @test all(isfinite, modal.Y.values[:, :, index])
            end
            rebuilt = compute(ModalTransformationProblem(modal))
            @test rebuilt.Z.values ≈ impedance rtol=1e-12
            @test rebuilt.Y.values ≈ admittance rtol=1e-12
        end
    end
end
