@testitem "LineModels / wideband and modal result records" tags=[:unit, :engine] begin
    import LineCableModels.LineModels as LM
    import LineCableModels.Units as U

    admittance = (poles=ComplexF64[-1, -2], residues=reshape(ComplexF64.(1:8), 2, 2, 2), direct=zeros(2, 2))
    groups = [(delay=1e-6, poles=ComplexF64[-3], residues=ones(ComplexF64, 2, 2, 1)),
        (delay=2e-6, poles=ComplexF64[-4, -5], residues=fill(2.0+0im, 2, 2, 2))]
    payload = WidebandModel(admittance, groups)
    errors = (rms=fill(0.001, 2, 2), max=fill(0.002, 2, 2), f_max=fill(100.0, 2, 2))
    passive = (checked=false, passive=nothing, bands=Tuple{Float64, Float64}[], worst=nothing,
        f_worst=nothing, perturbation=nothing)
    dc = (targets=nothing, achieved=nothing)
    retained = ComputationDetails(; formula=:ulm)
    result = LineModel(payload, ["west", "east"], 1000.0, [10.0, 100.0],
        (admittance=errors, propagation=deepcopy(errors)), passive, dc, retained)
    @test payload isa WidebandModel{Float64}
    @test poles(result) == [admittance.poles, groups[1].poles, groups[2].poles]
    @test residues(result) == [admittance.residues, groups[1].residues, groups[2].residues]
    @test delays(result) == [1e-6, 2e-6]
    @test fitting_error(result).admittance === errors
    @test passivity(result) === passive
    @test details(result) === retained
    @test computation_details(LineModelFormulation, result) === retained
    @test fieldnames(typeof(result)) == (:model, :coordinates, :line_length, :frequencies,
        :fitting_error, :passivity, :dc, :details)
    @test fieldnames(typeof(payload)) == (:admittance, :propagation)
    sections = [(constant=1.0, poles=[-1.0], residues=[2.0], delay=1e-6)]
    modal = ModalModel([1.0 0.0; 0.0 1.0], 50.0, sections, deepcopy(sections))
    @test modal isa ModalModel{Float64}
    @test poles(modal) == [[-1.0], [-1.0]]
    @test residues(modal) == [[2.0], [2.0]]
    @test delays(modal) == [1e-6]
    @test U.native_unit(U.Quantity{:pole}()) == U.units(:base, :dimensionless; per=(:base, :second))
    @test U.native_unit(U.Quantity{:residue}()) == U.units(:base, :dimensionless)
    @test U.native_unit(U.Quantity{:delay}()) == U.units(:base, :second)
    @test U.scale_factor(U.native_unit(U.Quantity{:fitting_error}()),
        U.display_unit(U.Quantity{:fitting_error}())) == 1.0
end

@testitem "LineModels / fitted evaluation at complex frequency and DC" tags=[:unit, :engine] begin
    for T in (Float64, Float32)
        s = Complex{T}[0, 3 + 4im, -1 + 2im, im]
        errors = (rms=zeros(T, 2, 2), max=zeros(T, 2, 2), f_max=zeros(T, 2, 2))
        passive = (checked=false, passive=nothing, bands=Tuple{T, T}[], worst=nothing,
            f_worst=nothing, perturbation=nothing)
        admittance = (poles=Complex{T}[-1, -2],
            residues=cat(Complex{T}[1 2; 3 4], Complex{T}[2 0; 1 3]; dims=3),
            direct=T[2 1; 0 3])
        groups = [(delay=T(0.1), poles=Complex{T}[-3],
                residues=reshape(Complex{T}[1 2; 3 4], 2, 2, 1)),
            (delay=T(0.2), poles=Complex{T}[-4],
                residues=reshape(Complex{T}[2 1; 0 3], 2, 2, 1))]
        wideband = LineModel(WidebandModel(admittance, groups), ["west", "east"],
            one(T), T[1, 2], (admittance=errors, propagation=errors), passive,
            (targets=nothing, achieved=nothing), ComputationDetails(; formula=:ulm))
        yc = @inferred Yc(wideband, s)
        h = @inferred H(wideband, s)
        @test eltype(yc) == eltype(h) == Complex{T}
        for (k, point) in enumerate(s)
            expected_yc = [2 + 1/(point+1) + 2/(point+2) 1 + 2/(point+1);
                3/(point+1) + 1/(point+2) 3 + 4/(point+1) + 3/(point+2)]
            expected_h = exp(-point*T(0.1)) / (point+3) * [1 2; 3 4] +
                exp(-point*T(0.2)) / (point+4) * [2 1; 0 3]
            @test yc[:, :, k] ≈ expected_yc
            @test h[:, :, k] ≈ expected_h
        end
        @test size(Yc(wideband, Complex{T}[])) == size(H(wideband, Complex{T}[])) == (2, 2, 0)

        # Nonorthogonal current modes distinguish transpose from inverse and voltage from current.
        impedance = [(constant=T(3), poles=T[-2], residues=T[12], delay=zero(T)),
            (constant=T(5), poles=T[-4], residues=T[20], delay=zero(T))]
        propagation = [(constant=zero(T), poles=T[-2, -4], residues=T[1, -1], delay=T(0.1)),
            (constant=T(0.7), poles=T[], residues=T[], delay=T(0.2))]
        modal = LineModel(ModalModel(T[1 2; 0 1], one(T), impedance, propagation),
            ["west", "east"], one(T), T[1, 2], (admittance=errors, propagation=errors),
            passive, (targets=nothing, achieved=nothing), ComputationDetails(; formula=:marti))
        yc = @inferred Yc(modal, s)
        h = @inferred H(modal, s)
        @test eltype(yc) == eltype(h) == Complex{T}
        for (k, point) in enumerate(s)
            y1, y2 = (point+2)/(3*(point+6)), (point+4)/(5*(point+8))
            h1 = exp(-point*T(0.1)) * 2/((point+2)*(point+4))
            h2 = T(0.7) * exp(-point*T(0.2))
            @test yc[:, :, k] ≈ [y1+4y2 2y2; 2y2 y2]
            @test h[:, :, k] ≈ [h1 2*(h2-h1); 0 h2]
        end
        @test size(Yc(modal, Complex{T}[])) == size(H(modal, Complex{T}[])) == (2, 2, 0)
    end
end
