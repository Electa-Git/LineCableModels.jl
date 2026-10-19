@testitem "ModalAnalysis / alignment / frequency-dependent character and phases" tags=[:unit,:engine] begin
    using LinearAlgebra: I
    import LineCableModels as LCM
    frequencies = [50.0,500.0,5000.0]
    current = repeat(reshape(Matrix{ComplexF64}(I,2,2),2,2,1),1,1,3)
    voltage = repeat(reshape(ComplexF64[2 .3im; .2 3],2,2,1),1,1,3)
    roots = ComplexF64[.01+.1im .02+.2im .03+.3im; .04+.4im .05+.5im .06+.6im]
    impedance = cat(ComplexF64[2+3im .4+.5im; .6+.7im 4+5im],
        ComplexF64[3+4im .5+.6im; .7+.8im 5+6im],
        ComplexF64[4+5im .6+.7im; .8+.9im 6+7im];dims=3)
    admittance = impedance .* (2+im)
    reference = LineParameters(ModalDomain(ModalOperators(voltage,current),roots),
        impedance,admittance,frequencies)
    order = [1 2 1; 2 1 2]
    phases = ComplexF64[1 im cis(.3); -im -1 cis(-.7)]
    trial_current, trial_voltage = similar(current), similar(voltage)
    trial_roots = similar(roots)
    trial_impedance, trial_admittance = similar(impedance), similar(admittance)
    for k in 1:3
        for mode in 1:2
            original = order[mode,k]
            trial_current[:,mode,k] = current[:,original,k] .* phases[mode,k]
            trial_voltage[:,mode,k] = voltage[:,original,k] .* phases[mode,k]
            trial_roots[mode,k] = roots[original,k]
            for row in 1:2
                ratio = phases[mode,k]/phases[row,k]
                trial_impedance[row,mode,k] = impedance[order[row,k],original,k]*ratio
                trial_admittance[row,mode,k] = admittance[order[row,k],original,k]*ratio
            end
        end
    end
    trial = LineParameters(ModalDomain(ModalOperators(trial_voltage,trial_current),trial_roots),
        trial_impedance,trial_admittance,frequencies)
    buffers = @inferred LCM.Commons.initialize_buffers(ModalAnalysis.align_modes,Float64,trial,(;),(;))
    aligned = @inferred ModalAnalysis.align_modes(trial,reference;buffers)
    independent = @inferred ModalAnalysis.align_modes(trial,reference)
    @test aligned.permutation == order
    @test aligned.overlap ≈ ones(2,3)
    @test aligned.permutation === buffers.permutation
    @test aligned.overlap === buffers.overlap
    @test Ti(aligned.parameters) ≈ current
    @test Tv(aligned.parameters) ≈ voltage
    @test gamma(aligned.parameters) == roots
    @test Z(aligned.parameters) ≈ impedance
    @test Y(aligned.parameters) ≈ admittance
    @test Z(transform(PhaseDomain,aligned.parameters)) ≈ Z(transform(PhaseDomain,reference))
    @test Y(transform(PhaseDomain,aligned.parameters)) ≈ Y(transform(PhaseDomain,reference))
    @test Ti(independent.parameters) !== buffers.Ti
    @test independent.permutation !== buffers.permutation
    @test validate(trial,ModalAnalysis.align_modes,reference) === trial
    ModalAnalysis.align_modes(reference,reference;buffers)
    @test independent.permutation == order
    @test aligned.permutation == repeat([1,2],1,3)
    mismatched = LineParameters(reference.domain,impedance,admittance,[60.0,500.0,5000.0])
    @test_throws DimensionMismatch ModalAnalysis.align_modes(trial,mismatched)
    damaged = deepcopy(trial)
    Ti(damaged)[:,1,2] .= 0
    @test_throws ArgumentError ModalAnalysis.align_modes(damaged,reference)
end

@testitem "ModalAnalysis / alignment / degenerate spans and unitary rotation" tags=[:unit,:engine] begin
    using LinearAlgebra: I
    tolerance = sqrt(eps(Float64))
    frequencies = [50.0,500.0]
    current = repeat(reshape(ComplexF64[1 .15 .2; 0 1 .1; 0 0 1],3,3,1),1,1,2)
    voltage = current .* (2+im)
    roots = ComplexF64[.1+.2im .1+.2im; .3+.4im .3+.4im; .1+.2im (.1+.2im)*sqrt(1+2tolerance)]
    impedance = repeat(reshape(ComplexF64[2+im .1 .2; .3 3+im .4; .5 .6 4+im],3,3,1),1,1,2)
    admittance = impedance .* (3-im)
    reference = LineParameters(ModalDomain(ModalOperators(voltage,current),roots),
        impedance,admittance,frequencies)
    order = [2,3,1]
    input_current,input_voltage = similar(current),similar(voltage)
    input_impedance,input_admittance = similar(impedance),similar(admittance)
    input_roots = 1.01 .* roots[order,:]
    for k in 1:2
        rotation = Matrix{ComplexF64}(I,3,3)
        rotation[2,2] = cis(.4)
        if k == 1
            rotation[[1,3],[1,3]] = [cos(.65) -sin(.65); sin(.65) cos(.65)]
        end
        input_current[:,:,k] = (current[:,:,k]*rotation)[:,order]
        input_voltage[:,:,k] = (voltage[:,:,k]*rotation)[:,order]
        input_impedance[:,:,k] = (rotation'*impedance[:,:,k]*rotation)[order,order]
        input_admittance[:,:,k] = (rotation'*admittance[:,:,k]*rotation)[order,order]
    end
    input = LineParameters(ModalDomain(ModalOperators(input_voltage,input_current),input_roots),
        input_impedance,input_admittance,frequencies)
    aligned = @inferred ModalAnalysis.align_modes(input,reference)
    @test aligned.cluster_sizes == [2 1; 1 1; 2 1]
    @test aligned.overlap ≈ ones(3,2)
    @test Ti(aligned.parameters) ≈ current
    @test Tv(aligned.parameters) ≈ voltage
    @test Z(aligned.parameters) ≈ impedance
    @test Y(aligned.parameters) ≈ admittance
    @test gamma(aligned.parameters) == hcat((input_roots[aligned.permutation[:,k],k] for k in 1:2)...)
    @test Z(transform(PhaseDomain,aligned.parameters)) ≈ Z(transform(PhaseDomain,input))
    @test Y(transform(PhaseDomain,aligned.parameters)) ≈ Y(transform(PhaseDomain,input))
    reference32 = LineParameters(ModalDomain(ModalOperators(ComplexF32.(voltage[:,:,1:1]),
        ComplexF32.(current[:,:,1:1])),ComplexF32.(roots[:,1:1])),
        ComplexF32.(impedance[:,:,1:1]),ComplexF32.(admittance[:,:,1:1]),Float32[50])
    input32 = LineParameters(ModalDomain(ModalOperators(ComplexF32.(input_voltage[:,:,1:1]),
        ComplexF32.(input_current[:,:,1:1])),ComplexF32.(input_roots[:,1:1])),
        ComplexF32.(input_impedance[:,:,1:1]),ComplexF32.(input_admittance[:,:,1:1]),Float32[50])
    aligned32 = @inferred ModalAnalysis.align_modes(input32,reference32)
    @test aligned32.cluster_sizes == reshape([2,1,2],3,1)
    @test eltype(Ti(aligned32.parameters)) === ComplexF32
    @test Ti(aligned32.parameters) ≈ Ti(reference32)
    @test Z(transform(PhaseDomain,aligned32.parameters)) ≈ Z(transform(PhaseDomain,input32))
    identity = reshape(Matrix{ComplexF64}(I,3,3),3,3,1)
    chain_roots = reshape(complex.(sqrt.([1.0,1+.75tolerance,1+1.5tolerance])),3,1)
    chain = LineParameters(ModalDomain(ModalOperators(identity,identity),chain_roots),identity,identity,[50.0])
    @test ModalAnalysis.align_modes(chain,chain).cluster_sizes == fill(3,3,1)
    zero_roots = reshape(ComplexF64[0,0,1],3,1)
    zeros_pair = LineParameters(ModalDomain(ModalOperators(identity,identity),zero_roots),identity,identity,[50.0])
    @test ModalAnalysis.align_modes(zeros_pair,zeros_pair).cluster_sizes == reshape([2,2,1],3,1)
end
