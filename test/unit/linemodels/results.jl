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
